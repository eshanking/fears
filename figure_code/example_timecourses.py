import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter, dir_manager
import pandas as pd

fig,ax = plt.subplots(nrows=3,ncols=2,figsize=(5,5))
labelsize = 10
#%% ROC data
suffix = '11012021_0000'
data_folder = 'results_' + suffix
exp_info_file = 'experiment_info_' + suffix + '.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
k_abs = exp_info.slopes

# exp_folders.reverse()
# k_abs = np.flip(k_abs)

exp_num = 1
pop_roc = exp_info.populations[0]

exp = exp_folders[exp_num]

# num_survived += 1
# num_extinct = 1

sim_files = os.listdir(path=exp)
sim_files = sorted(sim_files)

# find one extinct and one survived simulation

k = 0
found_extinct = False
found_survived = False

while found_extinct == False or found_survived == False:
    sim = sim_files[k]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    data = data[:,0:-1]
    data_t = data[-1,:]
    if any(data_t >= 1):
        num_survived = k
        found_survived = True
    else:
        num_extinct = k
        found_extinct = True
    k+=1

# plot survived timecourse

sim = sim_files[num_survived]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-1]
data = data[:,0:-1]
# data = data/np.max(data)
data_t = data[-1,:]
tcax = ax[0,0]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1,
               'linestyle':'-'
               # 'label':'Drug Concentration ($\u03BC$M)'
               }

data = data/np.max(data)

tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            data,
                                            tcax,
                                            drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=False,
                                            linewidth=1.5,
                                            drug_curve_label = '',
                                            labelsize = labelsize
                                            )
drug_ax.set_yticklabels('')
drug_ax.set_yticks([])
# plot extinct timecourse

sim = sim_files[num_extinct]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-1]
data = data[:,0:-1]
# data = data/np.max(data)
data_t = data[-1,:]
tcax = ax[0,1]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1
               # 'label':'Drug Concentration ($\u03BC$M)'
               }

data = data/np.max(data)

tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            data,
                                            tcax,
                                            drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=False,
                                            linewidth=1.5,
                                            labelsize = labelsize
                                            )

#%%

suffix = '11042021_0000'
data_folder = 'results_' + suffix
exp_info_file = 'experiment_info_' + suffix + '.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
# p_forget = exp_info.p_forget

# exp_folders.reverse()
# k_abs = np.flip(k_abs)

exp_num = 2
pop = exp_info.populations[0]
n_timestep = pop.n_timestep
exp = exp_folders[exp_num]

# num_survived += 1
# num_extinct = 1

sim_files = os.listdir(path=exp)
sim_files = sorted(sim_files)

# find one extinct and one survived simulation

k = 0
found_extinct = False
found_survived = False

while found_extinct == False or found_survived == False:
    sim = sim_files[k]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    data = data[:,0:-2]
    data_t = data[-1,:]
    if any(data_t >= 1):
        num_survived = k
        found_survived = True
    else:
        num_extinct = k
        found_extinct = True
    k+=1

# plot survived timecourse

sim = sim_files[num_survived]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-2]
survived_schedule = data[:,-1]

tcax = ax[1,0]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1
               # 'label':'Drug Concentration ($\u03BC$M)'
               }
tc = data[:,0:-2]
tc = tc/np.max(tc)
tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            tc,
                                            tcax,
                                            drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=True,
                                            legend_size=16,
                                            linewidth=1.5,
                                            drug_curve_label = '',
                                            labelsize = labelsize
                                            )
drug_ax.set_yticklabels('')
drug_ax.set_yticks([])
# plot extinct timecourse
#%%
sim = sim_files[num_extinct]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-2]
extinct_schedule = data[:,-1]

tcax = ax[1,1]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1
               # 'label':'Drug Concentration ($\u03BC$M)'
               }

tc = data[:,0:-2]
tc = tc/np.max(tc)

tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            tc,
                                            tcax,
                                            drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=False,
                                            linewidth=1.5,
                                            labelsize = labelsize,
                                            legend_size=16
                                            )

#%%

# plot dose times

x = np.arange(n_timestep-1)

timescale = pop.dose_schedule/pop.timestep_scale

# indx = np.argwhere(survived_schedule==1)
# indx = indx*timescale
# indx = indx[:,0].astype('int')
# survived_schedule = np.zeros(len(survived_schedule))
# survived_schedule[indx] = 1
#%%
# es_t = extinct_schedule
# indx = np.argwhere(extinct_schedule==1)
# indx = indx*timescale
# indx = indx[:,0].astype('int')
# extinct_schedule = np.zeros(len(extinct_schedule))
# extinct_schedule[indx] = 1

ax[2,0].plot(x,survived_schedule,linewidth=0.5,color='black')
ax[2,1].plot(x,extinct_schedule,linewidth=0.5,color='black')

#%% Adjust positions

pady = -0.08

ax[1,0] = plotter.shifty(ax[1,0],pady)
ax[1,1] = plotter.shifty(ax[1,1],pady)
ax[2,0] = plotter.shifty(ax[2,0],pady+0.01)
ax[2,1] = plotter.shifty(ax[2,1],pady+0.01)

padx = 0
ax[0,1] = plotter.shiftx(ax[0,1],padx)
ax[1,1] = plotter.shiftx(ax[1,1],padx)
ax[2,1] = plotter.shiftx(ax[2,1],padx)

# shrink regimen axes

shrink_pad = 0.2

p = ax[2,0].get_position()
p.y0 = p.y0+shrink_pad
ax[2,0].set_position(p)
p = ax[2,1].get_position()
p.y0 = p.y0+shrink_pad
ax[2,1].set_position(p)

#%% Adjust ticks and tick labels

# dose regimen plots
ax[2,0] = plotter.x_ticks_to_days(pop,ax[2,0])
ax[2,1] = plotter.x_ticks_to_days(pop,ax[2,1])

ax[1,0] = plotter.x_ticks_to_days(pop,ax[1,0])
ax[1,1] = plotter.x_ticks_to_days(pop,ax[1,1])

ax[0,0] = plotter.x_ticks_to_days(pop_roc,ax[0,0])
ax[0,1] = plotter.x_ticks_to_days(pop_roc,ax[0,1])
ax[0,0].set_xlabel('Days')
ax[0,1].set_xlabel('Days')

ax[1,0].set_xticklabels('')
ax[1,1].set_xticklabels('')
ax[2,1].set_yticks([])
ax[2,0].set_yticks([])
ax[2,1].spines['top'].set_visible(False)
ax[2,0].spines['top'].set_visible(False)

ax[1,1].set_yticklabels('')
ax[0,1].set_yticklabels('')
ax[1,1].set_yticks([])
ax[0,1].set_yticks([])

ax[2,0].set_xlabel('Days')
ax[2,1].set_xlabel('Days')

ax[0,0].set_ylabel('Proportion of max\n cell count')
ax[1,0].set_ylabel('Proportion of max\n cell count')

ax[1,0].legend(ncol=5,frameon=False,loc=(0,-1.3),fontsize=8)

results_manager.save_fig(fig,'example_timecourses.pdf',bbox_inches='tight')