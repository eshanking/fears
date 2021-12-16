import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter, dir_manager

def make_figure(roc_exp,adh_exp):
    
    fig,ax = plt.subplots(nrows=3,ncols=4,figsize=(6,4))
    labelsize = 10
    #%% ROC data
    # suffix = '11012021_0000' # lab machine
    # # suffix = '10272021_0001' # macbook
    # data_folder = 'results_' + suffix
    # exp_info_file = 'experiment_info_' + suffix + '.p'
    
    data_folder = roc_exp.results_path
    exp_info_file = roc_exp.experiment_info_path
    
    exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                                 exp_info_file)
    max_cells = exp_info.populations[0].max_cells
    n_sims = exp_info.n_sims
    k_abs = exp_info.slopes
    
    exp_num = 1
    pop_roc = exp_info.populations[0]
    
    exp = exp_folders[exp_num]
    
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
    # plot ROC survived timecourse
    
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
    drug_ax.plot(dc,color='black',linewidth=1)
    # drug_ax.set_yticklabels('')
    # drug_ax.set_yticks([])
    
    #%% plot ROC extinct timecourse
    
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
    drug_ax.plot(dc,color='black',linewidth=1)
    #%% nonadherance data
    
    # suffix = '11042021_0000' # lab machine
    # exp_num = 2 # lab machine
    # # exp_num = 0 # macbook
    # # suffix = '11082021_0000' # macbook
    # data_folder = 'results_' + suffix
    # exp_info_file = 'experiment_info_' + suffix + '.p'
    
    data_folder = adh_exp.results_path
    exp_info_file = adh_exp.experiment_info_path
    
    exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                                 exp_info_file)
    max_cells = exp_info.populations[0].max_cells
    n_sims = exp_info.n_sims
    
    pop = exp_info.populations[0]
    n_timestep = pop.n_timestep
    exp = exp_folders[exp_num]
    
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
    
    #%% plot nonadherance survived timecourse
    
    sim = sim_files[num_survived]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data[:,-2]
    survived_schedule = data[:,-1]
    
    tcax = ax[0,2]
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
    drug_ax = ax[1,2]
    drug_ax.plot(dc,color='black',linewidth=1)
    
    #%% plot nonadherance extinct timecourse
    
    sim = sim_files[num_extinct]
    sim = exp + os.sep + sim
    data = results_manager.get_data(sim)
    dc = data[:,-2]
    extinct_schedule = data[:,-1]
    
    tcax = ax[0,3]
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
    drug_ax = ax[1,3]
    drug_ax.plot(dc,color='black',linewidth=1)
    #%%  plot dose times
    
    x = np.arange(n_timestep-1)
    
    timescale = pop.dose_schedule/pop.timestep_scale
    
    ax[2,2].plot(x,survived_schedule,linewidth=0.2,color='black')
    ax[2,3].plot(x,extinct_schedule,linewidth=0.2,color='black')
    
    #%% Adjust positions
    
    for col in range(0,4):
        ax[0,col] = plotter.shrinky(ax[0,col],0.04)
        ax[1,col] = plotter.shrinky(ax[1,col],0.11)
        ax[1,col] = plotter.shifty(ax[1,col],0.04)
        ax[1,col].ticklabel_format(style='scientific',axis='y',
                                              scilimits=(0,3))
        
        ax[2,col] = plotter.shrinky(ax[2,col],0.2)
        ax[2,col] = plotter.shifty(ax[2,col],0.2)
    
    # shift right column to the right to make room for axis labels
    for col in range(2,4):
        for row in range(0,3):
            ax[row,col] = plotter.shiftx(ax[row,col],0.05)
            
    ax[2,2].spines['top'].set_visible(False)
    ax[2,3].spines['top'].set_visible(False)
    ax[2,2].spines['left'].set_visible(False)
    ax[2,3].spines['left'].set_visible(False)
    ax[2,2].spines['right'].set_visible(False)
    ax[2,3].spines['right'].set_visible(False)
    
    for row in range(0,3):
        ax[row,1].set_yticks([])
        ax[row,3].set_yticks([])
    
    ax[2,2].set_yticks([])
    
    ax[2,0].remove()
    ax[2,1].remove()
    
    #%% Adjust x axis tick labels
    
    xt = ax[0,0].get_xticks()
    xl = ax[0,0].get_xticklabels()
    xlim = ax[0,0].get_xlim()
    
    for col in range(0,2):
        ax[1,col] = plotter.x_ticks_to_days(pop_roc,ax[1,col])
        ax[1,col].set_xticks(xt)
        ax[1,col].set_xticklabels(xl)
        ax[1,col].set_xlim(xlim)
        ax[1,col].set_xlabel('Days')
        
    xt = ax[0,2].get_xticks()
    xl = ax[0,2].get_xticklabels()
    xlim = ax[0,2].get_xlim()
    
    for col in range(2,4):
        for row in range(1,3):
            ax[row,col] = plotter.x_ticks_to_days(pop_roc,ax[row,col])
            ax[row,col].set_xticks(xt)
            ax[row,col].set_xticklabels(xl)
            ax[row,col].set_xlim(xlim)
        ax[2,col].set_xlabel('Days')
    
    #%% Add axis labels
    
    ax[0,0].set_ylabel('Proportion',fontsize=8)
    ax[1,0].set_ylabel('Concentration ($\u03BC$M)',fontsize=8)
    
    ax[0,2].set_ylabel('Proportion',fontsize=8)
    ax[1,2].set_ylabel('Concentration ($\u03BC$M)',fontsize=8)
    
    #%% add panel labels
    
    alphabet = ['a','b','c','d','e','f','g','h','i','j']
    
    k = 0
    
    for row in range(2):
        for col in range(2):
            ax[row,col].annotate(alphabet[k],(0.92,1.07),xycoords='axes fraction')
            k+=1
            
    for row in range(0,2):
        for col in range(2,4):
            ax[row,col].annotate(alphabet[k],(0.92,1.1),xycoords='axes fraction')
            k+=1
    
    ax[2,2].annotate(alphabet[k],(0.92,1.4),xycoords='axes fraction')
    k+=1
    ax[2,3].annotate(alphabet[k],(0.92,1.4),xycoords='axes fraction')
    
    results_manager.save_fig(fig,'example_timecourses.pdf',bbox_inches='tight')
    
if __name__ == '__main__':
    make_figure()                            