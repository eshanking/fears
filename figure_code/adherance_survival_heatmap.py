import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
# from scipy import stats
from matplotlib.patches import Patch
import pickle
import sys

data_folder = 'results_06032021_0002'
exp_info_file = 'experiment_info_06032021_0002.p'

sys.path.append('/Users/eshanking/repos/')
from fears.classes import experiment_class_raw

def int_to_binary(num, pad=4):
    return bin(num)[2:].zfill(pad)

def chi_square(n,y,p=None):
    k = len(y)
    if p is None:
        p = np.ones(k)*1/k # uniform distribution
    if len(p) == 1:
        p = np.ones(k)*p
    V = 0
    for i in range(k):
        V += ((y[i]-n*p[i])**2)/(n*p[i])

    return V

def make_resultspath_absolute(filename):
    p = '..' + os.sep + 'results' + os.sep + filename
    return p

###############################################################################
# generate figure and axes

# exp_info = Experiment() # this is just a hack to suppress a warning
exp_info_path = make_resultspath_absolute(exp_info_file)
results_dir = make_resultspath_absolute(data_folder)

exp_info = pickle.load(open(exp_info_path,'rb')) # load experiment info

p1 = exp_info.populations[0]
max_cells = p1.max_cells

# results_dir = os.getcwd() + os.sep + data_folder # folder containing all of the results
experiment_folders = sorted(os.listdir(path=results_dir)) # each of these folders corresponds to a different k_abs
experiment_folders = [x for x in experiment_folders if x != '.DS_Store']

n_params = len(experiment_folders)

# Full page figure with room for caption
fig,ax = plt.subplots(2,1,figsize=(6.25,7.75),sharex=True)

# ax = fig.add_subplot()

##############################################################################
# data analysis and plotting

gap = int(p1.dose_schedule/p1.timestep_scale)
n_scheduled_doses = int(np.ceil(p1.n_timestep/gap))

exp_num = 2
exp = experiment_folders[exp_num]

p_drop = exp[exp.find('=')+1:]
p_drop = p_drop.replace(',', '.')

exp_path = results_dir + os.sep + exp
sim_files = os.listdir(path=exp_path)

n_sims = len(sim_files)

# resistant_regimen = np.zeros((1,n_scheduled_doses))
# extinct_regimen = np.zeros((1,n_scheduled_doses))
first_perished = True
first_survived = True

indx = 0

while indx < n_sims:
    
    sim = sim_files[indx]
    indx+=1
    
    sim_path = exp_path + os.sep + sim
    data_df = pd.read_csv(sim_path)
    data = data_df.to_numpy()
    counts = data[:,0:-3]
    regimen = data[0:n_scheduled_doses,-1]
    regimen = np.array([regimen])
    
    if any(counts[-1,:]>0.1*max_cells):
        if first_survived is True:
            # resistant_regimen[0,:] = regimen
            resistant_regimen = regimen
            first_survived = False
        else:
            resistant_regimen = np.concatenate((resistant_regimen,regimen),axis=0)
    else:
        if first_perished is True:
            # extinct_regimen[0,:] = regimen
            extinct_regimen = regimen
            first_perished = False
        else:
            extinct_regimen = np.concatenate((extinct_regimen,regimen),axis=0)   

# need to fix this
extinct_regimen = extinct_regimen[:,:-1]
resistant_regimen = resistant_regimen[:,:-1]

cmap = mpl.colors.ListedColormap(['cornflowerblue','w'])
# total_regimen = np.concatenate((resistant_regimen,extinct_regimen))
# ax[0].imshow(resistant_regimen,cmap=cmap,aspect=0.3)

n_patients = extinct_regimen.shape[0] + resistant_regimen.shape[0]
n_doses = extinct_regimen.shape[1]

aspect = n_doses/n_patients
ax[0].imshow(resistant_regimen,cmap=cmap,aspect=aspect)
ax[1].imshow(extinct_regimen,cmap=cmap,aspect=aspect)
# ax[1].imshow(extinct_regimen,cmap=cmap,aspect=0.3)

ax0_pos = ax[0].get_position()
ax1_pos = ax[1].get_position()
# ax[1].set_position([ax0_pos.x0,
#                     # ax1_pos.y0-ax1_pos.height*(ax0_pos.width/ax1_pos.width-1.3),
#                     ax1_pos.y0,
#                     ax0_pos.width,
#                     # ax1_pos.height])
#                     ax1_pos.height*ax0_pos.width/ax1_pos.width])  

# ax[0].set_position([ax0_pos.x0,
#                     ax1_pos.y0+0.05+ax1_pos.height,
#                     ax0_pos.width,
#                     ax0_pos.height])

new_height = ax1_pos.height*ax0_pos.width/ax1_pos.width

ax[1].set_position([ax0_pos.x0,
                    ax0_pos.y0-new_height-0.02,
                    ax0_pos.width,
                    new_height])
ax1_pos = ax[1].get_position()

# xticks = np.arange(0,n_doses,2)  
# ax[1].set_xticks(xticks)
# xlabels = [str(x) for x in np.arange(1,20,2)]
# ax[1].set_xticklabels(xlabels)
ax[1].set_xlabel('Dose number',fontsize=15)
ax[1].set_ylabel('Extinct',fontsize=15)
ax[0].set_ylabel('Resistant',fontsize=15)        

legend_elements = [Patch(facecolor='cornflowerblue',label='Missed dose'),Patch(facecolor='w',edgecolor='black',label='Scheduled dose')]
# ax[1].legend(handles=legend_elements,loc=(.15,-0.14),ncol=2,edgecolor='w') 
# ax[1].legend(handles=legend_elements,loc=(.15,ax1_pos.y1),ncol=2,edgecolor='w')   

title = '$p_{forget}$ = ' + p_drop

ax[0].set_title(title,fontsize=15)     

p_drop_t = p_drop.replace('.', ',')

plt.savefig('dose_regimen_heatmap_p=' + p_drop_t + '.pdf',bbox_inches="tight")
        
fig2, ax2 = plt.subplots(2,1,sharex=True,figsize=(3,4))

# for i in range(extinct_regimen.shape[0]):
#     y = np.correlate(extinct_regimen[i,:],extinct_regimen[i,:],mode='full')
#     y = y[int(np.floor(y.size/2)):]
#     x = np.arange(y.size)*2
#     ax2.scatter(x,y,color='red')
    
# for i in range(resistant_regimen.shape[0]):
#     y = np.correlate(resistant_regimen[i,:],resistant_regimen[i,:],mode='full')
#     y = y[int(np.floor(y.size/2)):]
#     x = np.arange(y.size)*2+0.5
#     ax2.scatter(x,y,color='blue')

# ax2.set_ylim(0.3,1)
# ax2.set_xlim(0,11)

n_survived = resistant_regimen.shape[0]
n_perished = extinct_regimen.shape[0]

survived_hist = np.sum(resistant_regimen,axis=0)
survived_hist = survived_hist/n_survived
perished_hist = np.sum(extinct_regimen,axis=0)
perished_hist = perished_hist/n_perished

dose_num = np.arange(len(survived_hist))
p = (1-float(p_drop))*np.ones(len(dose_num))

ax2[0].bar(dose_num,survived_hist,width=1,color='red',alpha=0.5,label='Resistant')
ax2[0].plot(dose_num,p,'--',color='black',label='$p_{remember}$ = 0.8')
ax2[1].bar(dose_num,perished_hist,width=1,color='blue',alpha=0.5,label='Extinct')
ax2[1].plot(dose_num,p,'--',color='black',label='$p_{remember}$ = 0.8')

# ax2[0].legend(loc=(0.25,0.05))
# ax2[1].legend(loc=(0.25,0.05))   

ax2[0].set_title('Resistant')
ax2[0].set_ylabel('Probability')
ax2[1].set_title('Extinct')
ax2[1].set_ylabel('Probability')    
ax2[1].set_xlabel('Dose number') 

plt.savefig('dose_regimen_histogram_p=' + p_drop_t + '.pdf',bbox_inches="tight")
# # p = p[0:10]
# # v_survived = chi_square(n_survived,np.sum(resistant_regimen[:,0:10],axis=0),p) 
# # v_perished = chi_square(n_perished,np.sum(extinct_regimen[:,0:10],axis=0),p)
# p_value = []

# n_taken_survived = np.sum(resistant_regimen,axis=0)
# n_taken_perished = np.sum(extinct_regimen,axis=0)
# n_missed_survived = n_survived - n_taken_survived
# n_missed_perished = n_perished - n_taken_perished

# for n in dose_num:
#     table = [[n_taken_survived[n],n_missed_survived[n]],[n_taken_perished[n],n_missed_perished[n]]]
#     oddsratio, p_value_t = stats.fisher_exact(table)
#     p_value.append(p_value_t)
    