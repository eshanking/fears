import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter, dir_manager

fig,ax = plt.subplots(nrows=5,ncols=2,figsize=(4,5))
labelsize = 10
#%% ROC data
# suffix = '11012021_0000' # lab machine
suffix = '10272021_0001' # macbook
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
#%%
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
label_kwargs = {'align':False}
select_labels = [6,14,15]
label_xpos = [1200,2500,4500]

data = data/np.max(data)

tcax,t = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            data,
                                            tcax,
                                            # drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=True,
                                            linewidth=1.5,
                                            drug_curve_label = '',
                                            labelsize = labelsize,
                                            label_lines = True,
                                            select_labels = select_labels,
                                            label_xpos=label_xpos,
                                            label_kwargs=label_kwargs
                                            )
drug_ax = ax[1,0]
drug_ax.plot(dc,color='black')
# drug_ax.set_yticklabels('')
# drug_ax.set_yticks([])

#%% plot extinct timecourse

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

tcax,t = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            data,
                                            tcax,
                                            # drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=False,
                                            linewidth=1.5,
                                            labelsize = labelsize
                                            )
drug_ax = ax[1,1]
drug_ax.plot(dc,color='black')
#%% nonadherance data

# suffix = '11042021_0000' # lab machine
# exp_num = 2 # lab machine
exp_num = 0 # macbook
suffix = '11082021_0000' # macbook
data_folder = 'results_' + suffix
exp_info_file = 'experiment_info_' + suffix + '.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
# p_forget = exp_info.p_forget

# exp_folders.reverse()
# k_abs = np.flip(k_abs)


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

#%% plot survived timecourse

sim = sim_files[num_survived]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-2]
survived_schedule = data[:,-1]

tcax = ax[2,0]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1
               # 'label':'Drug Concentration ($\u03BC$M)'
               }
label_kwargs = {'align':False}
tc = data[:,0:-2]
tc = tc/np.max(tc)

select_labels = [2,6]
label_xpos = [350,1100]

tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            tc,
                                            tcax,
                                            # drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=True,
                                            legend_size=16,
                                            linewidth=1.5,
                                            drug_curve_label = '',
                                            labelsize = labelsize,
                                            label_lines=True,
                                            select_labels=select_labels,
                                            label_xpos=label_xpos,
                                            label_kwargs=label_kwargs
                                            )
drug_ax = ax[3,0]
drug_ax.plot(dc,color='black')
# drug_ax.set_yticklabels('')
# drug_ax.set_yticks([])

#%% plot extinct timecourse

sim = sim_files[num_extinct]
sim = exp + os.sep + sim
data = results_manager.get_data(sim)
dc = data[:,-2]
extinct_schedule = data[:,-1]

tcax = ax[2,1]
drug_kwargs = {'alpha':1,
               'color':'black',
               'linewidth':1
               # 'label':'Drug Concentration ($\u03BC$M)'
               }

tc = data[:,0:-2]
tc = tc/np.max(tc)

select_labels = [0]
label_xpos = [100]

tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[exp_num],
                                            tc,
                                            tcax,
                                            # drug_curve=dc,
                                            drug_ax_sci_notation=True,
                                            drug_kwargs=drug_kwargs,
                                            legend_labels=True,
                                            linewidth=1.5,
                                            labelsize = labelsize,
                                            legend_size=16,
                                            select_labels=select_labels,
                                            label_xpos=label_xpos,
                                            label_lines=True
                                            )
drug_ax = ax[3,1]
drug_ax.plot(dc,color='black')
#%%  plot dose times

x = np.arange(n_timestep-1)

timescale = pop.dose_schedule/pop.timestep_scale

#%%

ax[4,0].plot(x,survived_schedule,linewidth=0.5,color='black')
ax[4,1].plot(x,extinct_schedule,linewidth=0.5,color='black')

#%% Adjust positions

ax[1,0] = plotter.shrinky(ax[1,0],0.05)
ax[1,1] = plotter.shrinky(ax[1,1],0.05)

ax[1,0] = plotter.shifty(ax[1,0],-0.01)
ax[1,1] = plotter.shifty(ax[1,1],-0.01)

ax[2,0] = plotter.shifty(ax[2,0],-0.08)
ax[2,1] = plotter.shifty(ax[2,1],-0.08)

ax[3,0] = plotter.shrinky(ax[3,0],0.04)
ax[3,1] = plotter.shrinky(ax[3,1],0.04)

ax[3,0] = plotter.shifty(ax[3,0],-0.11)
ax[3,1] = plotter.shifty(ax[3,1],-0.11)

ax[4,0] = plotter.shrinky(ax[4,0],0.1)
ax[4,1] = plotter.shrinky(ax[4,1],0.1)

ax[4,0] = plotter.shifty(ax[4,0],-0.08)
ax[4,1] = plotter.shifty(ax[4,1],-0.08)

#%% Adjust ticks and tick labels

# dose regimen plots

ax[0,0] = plotter.x_ticks_to_days(pop_roc,ax[0,0])
ax[0,1] = plotter.x_ticks_to_days(pop_roc,ax[0,1])

ax[1,0] = plotter.x_ticks_to_days(pop_roc,ax[1,0])
ax[1,1] = plotter.x_ticks_to_days(pop_roc,ax[1,1])

ax[1,0].ticklabel_format(style='scientific',axis='y',
                                     scilimits=(0,3))
ax[1,1].ticklabel_format(style='scientific',axis='y',
                                     scilimits=(0,3))

ax[2,0] = plotter.x_ticks_to_days(pop,ax[2,0])
ax[2,1] = plotter.x_ticks_to_days(pop,ax[2,1])

ax[3,0] = plotter.x_ticks_to_days(pop,ax[3,0])
ax[3,1] = plotter.x_ticks_to_days(pop,ax[3,1])

ax[3,0].ticklabel_format(style='scientific',axis='y',
                                     scilimits=(0,3))
ax[3,1].ticklabel_format(style='scientific',axis='y',
                                     scilimits=(0,3))

ax[4,0] = plotter.x_ticks_to_days(pop,ax[4,0])
ax[4,1] = plotter.x_ticks_to_days(pop,ax[4,1])

ax[4,1].spines['top'].set_visible(False)
ax[4,0].spines['top'].set_visible(False)

ax[1,1].set_yticklabels('')
ax[0,1].set_yticklabels('')
ax[1,1].set_yticks([])
ax[0,1].set_yticks([])
ax[2,1].set_yticklabels('')
ax[2,1].set_yticks([])
ax[3,1].set_yticklabels('')
ax[3,1].set_yticks([])
ax[4,1].set_yticklabels('')
ax[4,1].set_yticks([])
ax[4,0].set_yticklabels('')
ax[4,0].set_yticks([])

ax[1,0].set_xlabel('Days')
ax[1,1].set_xlabel('Days')

ax[4,0].set_xlabel('Days')
ax[4,1].set_xlabel('Days')

xt = ax[1,0].get_xticks()
xtl = ax[1,0].get_xticklabels()
xl = ax[1,0].get_xlim()

ax[0,0].set_xticks(xt)
ax[0,1].set_xticks(xt)
ax[0,0].set_xticklabels(xtl)
ax[0,1].set_xticklabels(xtl)
ax[0,0].set_xlim(xl)
ax[0,1].set_xlim(xl)

xt = ax[2,0].get_xticks()
xtl = ax[2,0].get_xticklabels()
xl = ax[2,0].get_xlim()

ax[3,0].set_xticks(xt)
ax[3,1].set_xticks(xt)
ax[3,0].set_xticklabels(xtl)
ax[3,1].set_xticklabels(xtl)
ax[3,0].set_xlim(xl)
ax[3,1].set_xlim(xl)

ax[4,0].set_xticks(xt)
ax[4,1].set_xticks(xt)
ax[4,0].set_xticklabels(xtl)
ax[4,1].set_xticklabels(xtl)
ax[4,0].set_xlim(xl)
ax[4,1].set_xlim(xl)

fontsize = 8
ax[0,0].set_ylabel('Proportion',fontsize=fontsize)
ax[2,0].set_ylabel('Proportion',fontsize=fontsize)

ax[1,0].set_ylabel('Drug \n concentration',fontsize=fontsize)
ax[3,0].set_ylabel('Drug \n concentration',fontsize=fontsize)

# ax[2,0].legend(ncol=4,frameon=False,loc=(0,-4.3),fontsize=8)

ax[0,0].set_title('Resistant')
ax[0,1].set_title('Extinct')

#%% Add caption annotations

cords = (0.95,1)
cords2 = (0.95,1.1)
cords3 = (0.95,1.35)

ax[0,0].annotate('a',xy = cords,xycoords='axes fraction')
ax[0,1].annotate('b',xy = cords,xycoords='axes fraction')
ax[1,0].annotate('c',xy = cords2,xycoords='axes fraction')
ax[1,1].annotate('d',xy = cords2,xycoords='axes fraction')
ax[2,0].annotate('e',xy = cords,xycoords='axes fraction')
ax[2,1].annotate('f',xy = cords,xycoords='axes fraction')
ax[3,0].annotate('g',xy = cords2,xycoords='axes fraction')
ax[3,1].annotate('h',xy = cords2,xycoords='axes fraction')
ax[4,0].annotate('i',xy = cords3,xycoords='axes fraction')
ax[4,1].annotate('j',xy = cords3,xycoords='axes fraction')

results_manager.save_fig(fig,'example_timecourses.pdf',bbox_inches='tight')