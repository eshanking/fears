from fears.utils import results_manager,plotter,dir_manager
import matplotlib.pyplot as plt
import numpy as np
import os

suffix = '07212021_0001'
data_folder = 'results_' + suffix
exp_info_file = 'experiment_info_' + suffix + '.p'

fig = plt.figure(figsize=(10,2),constrained_layout=True)
width_ratios = [1.5,1.5,1,1,1,1,1,1]
gs = fig.add_gridspec(nrows=2,ncols=8,width_ratios=width_ratios)
tc_axes = []
for i in range(2):
    for j in range(2):
        a = fig.add_subplot(gs[i,j])
        tc_axes.append(a)

km_axes = []        
km_axes.append( fig.add_subplot(gs[0:2,2:4]) )
km_axes.append( fig.add_subplot(gs[0:2,4:6]) )
km_axes.append( fig.add_subplot(gs[0:2,6:8]) )

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
k_abs = exp_info.slopes


# make a dummy barchart
# k_abs = k_abs[2:]
x = np.arange(len(k_abs))
barchart_data = np.ones(len(k_abs))*100

drug_axes = []

exp_folders.reverse()
k_abs = np.flip(k_abs)
pop = exp_info.populations[0]
for exp in exp_folders:
    k_abs_t = exp[exp.find('=')+1:]
    k_abs_t = float(k_abs_t)
    num = np.argwhere(k_abs == k_abs_t)
    num = num[0,0]
    
    # generate timecourse axes
    

    # da = tcax.twinx()
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    survive_count = 0
    counts_total = None
    
    death_event_obs = np.zeros(n_sims)
    death_event_times = np.zeros(n_sims)
    
    gen14_resistance_obs = np.zeros(n_sims)
    gen14_resistance_times = np.zeros(n_sims)
    
    gen15_resistance_obs = np.zeros(n_sims)
    gen15_resistance_times = np.zeros(n_sims)
    
    k=0
    # while k < len(sim_files):
    while k < 10:
    # for sim in sim_files:
        sim = sim_files[k]
        sim = exp + os.sep + sim
        data = results_manager.get_data(sim)
        dc = data[:,-1]
        data = data[:,0:-1]
        # data = data/np.max(data)
        data_t = data[-1,:]
        
        # check to see if any genotypes are at least 10% of the max cell count
        if any(data_t >= 1):
            survive_count += 1
            if counts_total is None:
                counts_total = data
            else:
                counts_total += data
        # data = data/np.max(data)
        
        # exp_info.populations[num].counts_log_scale = True
        
        
        death_event_obs[k],death_event_times[k] = \
            exp_info.extinction_time(pop,data,thresh=1)
            
        gen15_resistance_obs[k],gen15_resistance_times[k] = \
            exp_info.resistance_time(pop,data,15,thresh=.1)

        gen14_resistance_obs[k],gen14_resistance_times[k] = \
            exp_info.resistance_time(pop,data,14,thresh=.1)
        
        data = data/max_cells
        if k==0:
            drug_kwargs = {'alpha':0.7,
                           'color':'black',
                           'linewidth':2,
                           'label':'Drug Concentration ($\u03BC$M)'
                           }
            tc_axes[num],drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tc_axes[num],
                                                        drug_curve=dc,
                                                        drug_ax_sci_notation=True,
                                                        drug_kwargs=drug_kwargs,
                                                        legend_labels=False,
                                                        grayscale=True,
                                                        color='gray', 
                                                        linewidth=1,
                                                        labelsize=10,
                                                        alpha=0.7
                                                        )
            drug_ax.set_ylabel('')
            drug_axes.append( drug_ax )
        else:
            tc_axes[num],da = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tc_axes[num],
                                                        grayscale=True,
                                                        color='gray',
                                                        legend_labels=False,
                                                        linewidth=2,
                                                        labelsize=10,
                                                        alpha=0.2
                                                        )            
        # drug_ax.set_ylim(0,10**4)
        k+=1
        
    if survive_count > 0:
        counts_avg = counts_total/survive_count
        # counts_avg = counts_avg/np.max(counts_avg) 
        # counts_avg = counts_total
        counts_avg = counts_avg/np.max(counts_avg)
        tc_axes[num],temp = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                               counts_avg,
                                               tc_axes[num],
                                               labelsize=10)
    km_axes[0] = plotter.plot_kaplan_meier(pop,
                                      death_event_times,
                                      ax=km_axes[0],
                                      n_sims=n_sims,
                                      label=str(k_abs_t),
                                      mode='survival')
    
    km_axes[1] = plotter.plot_kaplan_meier(pop,
                                      gen14_resistance_times,
                                      ax=km_axes[1],
                                      n_sims=n_sims,
                                      label=str(k_abs_t),
                                      mode='resistant')
    
    km_axes[2] = plotter.plot_kaplan_meier(pop,
                                      gen15_resistance_times,
                                      ax=km_axes[2],
                                      n_sims=n_sims,
                                      label=str(k_abs_t),
                                      mode='resistant')
    
for da in drug_axes:
    da.ticklabel_format(style='sci',axis='y',scilimits=(0,4))
    
for ax in tc_axes:
    
    ax.set_box_aspect(1)
    # p = ax.get_position()
    # ax.set_position([p.x0,p.y0,p.width*5,p.height*5])
    
    

for ax in km_axes:
    ax.set_ylim([-5,105])