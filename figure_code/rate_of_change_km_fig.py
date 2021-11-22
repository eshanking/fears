import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter, dir_manager
import pandas as pd

def make_fig(roc_exp):
        # suffix = '11012021_0000' # lab machine
    # suffix = '10272021_0001' # macbook
    # data_folder = 'results_' + suffix
    # exp_info_file = 'experiment_info_' + suffix + '.p'
    
    data_folder = roc_exp.results_path
    exp_info_file = roc_exp.experiment_info_path
    
    exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                                 exp_info_file)
    max_cells = exp_info.populations[0].max_cells
    n_sims = exp_info.n_sims
    k_abs = exp_info.slopes
    
    # exp_folders.reverse()
    # k_abs = np.flip(k_abs)
    
    fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(8,2.5))
    
    pop = exp_info.populations[0]
    
    km_data = {'survival':{},
               'resistance 1110':{},
               'resistance 1111':{}}
    
    for exp in exp_folders:
    
        k_abs_t = exp[exp.find('=')+1:]
        k_abs_t = k_abs_t.replace(',','.')
        k_abs_t = float(k_abs_t)
        
        num = np.argwhere(k_abs == k_abs_t)
        num = num[0,0]
        
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        
        # KM data 
        death_event_obs = np.zeros(n_sims)
        death_event_times = np.zeros(n_sims)
        
        # time to genotype 2
        gen14_resistance_obs = np.zeros(n_sims)
        gen14_resistance_times = np.zeros(n_sims)
        
        # time to genotype 6
        gen15_resistance_obs = np.zeros(n_sims)
        gen15_resistance_times = np.zeros(n_sims)
        
        k=0
        while k < len(sim_files):
        # while k < 10:
            sim = sim_files[k]
            sim = exp + os.sep + sim
            data = results_manager.get_data(sim)
            dc = data[:,-1]
            data = data[:,0:-1]
            
            death_event_obs[k],death_event_times[k] = \
                exp_info.extinction_time(pop,data,thresh=1)
                
            gen14_resistance_obs[k],gen14_resistance_times[k] = \
                exp_info.resistance_time(pop,data,14,thresh=.1)
    
            gen15_resistance_obs[k],gen15_resistance_times[k] = \
                exp_info.resistance_time(pop,data,15,thresh=.1)
                
            k+=1
            
        ax[0] = plotter.plot_kaplan_meier(pop,
                                          death_event_times,
                                          ax=ax[0],
                                          n_sims=n_sims,
                                          label=str(k_abs_t),
                                          mode='survival')
        
        ax[1] = plotter.plot_kaplan_meier(pop,
                                          gen14_resistance_times,
                                          ax=ax[1],
                                          n_sims=n_sims,
                                          label=str(k_abs_t),
                                          mode='resistant')
        
        ax[2] = plotter.plot_kaplan_meier(pop,
                                          gen15_resistance_times,
                                          ax=ax[2],
                                          n_sims=n_sims,
                                          label=str(k_abs_t),
                                          mode='resistant')
        
        km_data['survival'][str(k_abs_t)] = death_event_times
        km_data['resistance 1110'][str(k_abs_t)] = gen14_resistance_times
        km_data['resistance 1111'][str(k_abs_t)] = gen15_resistance_times
    
    for a in ax:
        a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(False)
        
    ax[2].legend(frameon=False,loc='upper right',title='$k_{abs}$',fontsize=8)
    
    pad = 0.05
    
    ax[0].set_ylabel('% surviving')
    pos1 = ax[1].get_position()
    pos1.x0 = pos1.x0 + pad
    pos1.x1 = pos1.x1 + pad
    ax[1].set_position(pos1)
    
    pos2 = ax[2].get_position()
    pos2.x0 = pos2.x0 + pad*2
    pos2.x1 = pos2.x1 + pad*2
    ax[2].set_position(pos2)
    
    ax[0].set_title('Survival of infectious agent',fontsize=8)
    ax[1].set_title('Resistant genotype = 1110',fontsize=8)
    ax[2].set_title('Resistant genotype = 1111',fontsize=8)
    
    # max sure all x lims are the same
    xmax = ax[0].get_xlim()[1]
    for a in ax:
        if a.get_xlim()[1] > xmax:
            xmax = a.get_xlim()[1]
    
    for a in ax:
        a.set_xlim([0,xmax])
    #%%
    # fix ax[0] xlabels ¯\_(ツ)_/¯
    xt = ax[1].get_xticks()
    xl = ax[1].get_xticklabels()
    ax[0].set_xticks(xt)
    ax[0].set_xticklabels(xl)
    ax[0].set_xlim([0,xmax])
    
    # fix ax[2] xlabels ¯\_(ツ)_/¯
    ax[2].set_xticks(xt)
    ax[2].set_xticklabels(xl)
    ax[2].set_xlim([0,xmax])
    results_manager.save_fig(fig,'roc_km_curve.pdf',bbox_inches='tight')
    #%%
    # perform pairwise log-rank tests and compute p values
    analysis_keys = list(km_data.keys()) # endpoints being analyzed
    experiment_keys = [str(p) for p in k_abs] # experiments performed
    
    comparisons = [] # vector of all pairwise comparisons without duplicates
    
    for i in range(len(experiment_keys)):
        j = i+1
        while j < len(experiment_keys):
            pair = (k_abs[i],k_abs[j])
            j+=1
            comparisons.append(pair)
    
    p_values = {'survival':{},
                'resistance 1110':{},
                'resistance 1111':{}}
    
    n_tests = len(k_abs)-1
    
    for ak in  analysis_keys:
        for pair in comparisons:
            key0 = str(pair[0])
            key1 = str(pair[1])
            sr = exp_info.log_rank_test(km_data[ak][key0],km_data[ak][key1])
            p_values[ak][str(pair)] = float(sr.p_value)*n_tests # Mutliple hypothesis testing correction
    
    p_values = pd.DataFrame(p_values)
    result_path = dir_manager.make_resultspath_absolute(
        'rate_of_change_km_curves_p_values.csv')
    
    p_values.to_csv(result_path)
#%%
if __name__ == '__main__':
    make_figure()
#%%

def mutation_burden(pop,counts):
    
    return


          