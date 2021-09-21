import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter, dir_manager
import os

suffix = '07212021_0001'
data_folder = 'results_' + suffix
exp_info_file = 'experiment_info_' + suffix + '.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
k_abs = exp_info.slopes

exp_folders.reverse()
k_abs = np.flip(k_abs)

fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(4,4))

pop = exp_info.populations[0]

ax = ax.reshape((len(k_abs),))

axnum = 0
tc_axes=[]
drug_axes=[]

for exp in exp_folders:
    k_abs_t = exp[exp.find('=')+1:]
    k_abs_t = float(k_abs_t)
    num = np.argwhere(k_abs == k_abs_t)
    num = num[0,0]
    
    # generate timecourse axes
    
    tcax = ax[axnum]
    # da = tcax.twinx()
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    survive_count = 0
    counts_total = None
    
    k=0
    while k < len(sim_files):
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
        data = data/max_cells
        if k==0:
            drug_kwargs = {'alpha':0.7,
                           'color':'black',
                           'linewidth':2,
                           'label':'Drug Concentration ($\u03BC$M)'
                           }
            tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tcax,
                                                        drug_curve=dc,
                                                        drug_ax_sci_notation=True,
                                                        drug_kwargs=drug_kwargs,
                                                        legend_labels=False,
                                                        grayscale=True,
                                                        color='gray', 
                                                        linewidth=1,
                                                        labelsize=12,
                                                        alpha=0.7
                                                        )
            drug_ax.set_ylabel('')
            drug_axes.append( drug_ax )
        else:
            tcax,da = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tcax,
                                                        grayscale=True,
                                                        color='gray',
                                                        legend_labels=False,
                                                        linewidth=2,
                                                        labelsize=12,
                                                        alpha=0.2
                                                        )            
        # drug_ax.set_ylim(0,10**4)
        k+=1
        
    if survive_count > 0:
        counts_avg = counts_total/survive_count
        # counts_avg = counts_avg/np.max(counts_avg) 
        # counts_avg = counts_total
        counts_avg = counts_avg/np.max(counts_avg)
        tcax,temp = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                               counts_avg,
                                               tcax,
                                               labelsize=12)
    
    # t = np.arange(len(dc))
    # t = t*exp_info.populations[0].timestep_scale/24
    # da.plot(t,dc)
    
    tc_axes.append( tcax )
    axnum+=1