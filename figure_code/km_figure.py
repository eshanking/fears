import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter

data_folder = 'results_07202021_0000'
exp_info_file = 'experiment_info_07202021_0000.p'

fig,ax = plt.subplots(figsize=(4,3.5))
labelsize=12

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
p_drop = exp_info.prob_drops


# make a dummy barchart
# p_drop = p_drop[2:]
x = np.arange(len(p_drop))
barchart_data = np.ones(len(p_drop))*100
rects = ax.barh(x,barchart_data,color='slategrey',facecolor='w')

tc_axes = []
drug_axes = []
pop_axes = []

exp_folders.reverse()
p_drop = np.flip(p_drop)

thresh = 1
pop = exp_info.populations[0]

# data_extinct = np.zeros((999,1))
for exp in exp_folders:
    p_drop_t = exp[exp.find('=')+1:]
    p_drop_t = p_drop_t.replace(',','.')
    p_drop_t = float(p_drop_t)
    num = np.argwhere(p_drop == p_drop_t)
    num = num[0,0]
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    
    event_obs = np.zeros(n_sims)
    event_times = np.zeros(n_sims)
    
    k=0
    while k < len(sim_files):
    # while k < 10:
        sim = sim_files[k]
        sim = exp + os.sep + sim
        data = results_manager.get_data(sim)
        dc = data[:,-2]
        data = data[:,0:-2]
        event_obs[k],event_times[k] = \
            exp_info.extinction_time(pop,data,thresh=1)
        k+=1
        
    ax = plotter.plot_kaplan_meier(pop,event_times,ax=ax,n_sims=n_sims,
                                   label=str(p_drop_t),
                                   mode='extinct')
    
ax.legend(frameon=False,loc='upper left',title='$p_{forget}$')

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

results_manager.save_fig(fig,'km_curve.pdf',bbox_inches='tight')