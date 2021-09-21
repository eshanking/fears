import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import sys
import os

sys.path.append('/Users/eshanking/repos/')
from fears.classes import experiment_class_raw
from fears.src import utils

##############################################################################
data_folder = 'results_06032021_0002'
exp_info_file = 'experiment_info_06032021_0002.p'

exp_num = 2
##############################################################################
def int_to_binary(num, pad=4):
    return bin(num)[2:].zfill(pad)

def max_time_between_doses(regimen):
    last_dose = 0
    max_interval = 0
    for i in range(len(regimen)):
        if regimen[i] == 1:
            interval = i-last_dose
            if interval > max_interval:
                max_interval = interval
            last_dose = i
    return max_interval

##############################################################################
# generate figures and axes
fig,ax = plt.subplots()

# load data and metadata

exp_info_path = utils.make_resultspath_absolute(exp_info_file)
results_dir = utils.make_resultspath_absolute(data_folder)

exp_info = pickle.load(open(exp_info_path,'rb')) # load experiment info

p1 = exp_info.populations[0]

experiment_folders = sorted(os.listdir(path=results_dir)) # each of these folders corresponds to a different k_abs
experiment_folders = [x for x in experiment_folders if x != '.DS_Store'] # remove mac hidden folders

n_params = len(experiment_folders)

# compute results

gap = int(p1.dose_schedule/p1.timestep_scale)
n_scheduled_doses = int(np.ceil(p1.n_timestep/gap))
exp = experiment_folders[exp_num]

max_cells = p1.max_cells

p_drop = exp[exp.find('=')+1:]
p_drop = p_drop.replace(',', '.')

exp_path = results_dir + os.sep + exp
sim_files = os.listdir(path=exp_path)

n_sims = len(sim_files)

indx = 0
resistant_intervals = []
extinct_intervals = []

while indx < n_sims:
    sim = sim_files[indx]
    indx+=1
    
    sim_path = exp_path + os.sep + sim
    data_df = pd.read_csv(sim_path)
    data = data_df.to_numpy()
    counts = data[:,0:-3]
    regimen = data[0:n_scheduled_doses,-1]
    regimen = np.array([regimen])
    
    m = max_time_between_doses(regimen[0,:])
    
    if any(counts[-1,:]>0.1*max_cells):
        resistant_intervals.append(m)
    else:
        extinct_intervals.append(m)

binwidth = 2
# r_bins=range(min(resistant_intervals), max(resistant_intervals) + binwidth, binwidth)
# e_bins=range(min(extinct_intervals), max(extinct_intervals) + binwidth, binwidth)

r_bins=range(5, 40 + binwidth, binwidth)
e_bins = r_bins
# e_bins=range(min(extinct_intervals), max(extinct_intervals) + binwidth, binwidth)

ax.hist(resistant_intervals,alpha=0.5,density=True,bins=r_bins)
ax.hist(extinct_intervals,alpha=0.5,density=True,bins=e_bins)

ax.legend(['Resistant','Extinct'])
ax.set_ylabel('Probability',fontsize=15)
ax.set_xlabel('Max sequential missed doses',fontsize=15)
ax.set_title('$p_{forget}$ = ' + str(p_drop),fontsize=15)

plt.savefig('figures' + os.sep + 'max_interval_hist_p=' + p_drop + '.pdf',bbox_inches="tight")