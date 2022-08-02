import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager
import pickle

def km_curve(exp=None,exp_info_path=None,resistance_outcome=[14,15]):
    """Returns a dictionary of dictionaries of K-M curves from the 
       given experiment. Each experimental condition has two 
       resistance curves (defined by resistance_outcome) and a
       survival curve.

    Args:
        exp (fears Experiment object, optional): Experiment to analyze. Defaults to None.
        exp_info_path (str, optional): Optional path to Experiment object. Defaults to None.
        resistance_outcome (list, optional): List of resistance outcomes. Defaults to [14,15].

    Returns:
        dict: Dict of dicts containing KM curves. Sub-dictionaries are different
        experimental conditions.
    """

    if exp is None:
        exp =  pickle.load(open(exp_info_path,'rb'))

    exp_folders,exp_info = results_manager.get_experiment_results(exp=exp)

    n_sims = exp_info.n_sims
    
    pop = exp_info.populations[0]
    
    km_data = {}
    
    for exp in exp_folders:
    
        k_abs_t = exp[exp.find('=')+1:]
        k_abs_t = k_abs_t.replace(',','.')
        k_abs_t = float(k_abs_t)
        
        k_abs_t = round(k_abs_t,10)
        # print(f"{k_abs_t:.2e}")
        
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        
        # KM data 
        death_event_obs = np.zeros(n_sims)
        death_event_times = np.zeros(n_sims)
        
        # time to genotype 2
        # gen14_resistance_obs = np.zeros(n_sims)
        # gen14_resistance_times = np.zeros(n_sims)
        gen1_resistance_obs = np.zeros(n_sims)
        gen1_resistance_times = np.zeros(n_sims)
        
        # time to genotype 6
        gen2_resistance_obs = np.zeros(n_sims)
        gen2_resistance_times = np.zeros(n_sims)
        
        k=0

        km_data_t = {}
        while k < len(sim_files):

            sim = sim_files[k]
            sim = exp + os.sep + sim
            data_dict = results_manager.get_data(sim)

            data = data_dict['counts']
            
            death_event_obs[k],death_event_times[k] = \
                exp_info.extinction_time(pop,data,thresh=1)


            gen1_resistance_obs[k],gen1_resistance_times[k] = \
                exp_info.resistance_time(pop,data,resistance_outcome[0],thresh=.1)
    
            gen2_resistance_obs[k],gen2_resistance_times[k] = \
                exp_info.resistance_time(pop,data,resistance_outcome[1],thresh=.1)
                
            k+=1
            
        km_data_t['survival'] = death_event_times
        key1 = 'resistance ' + pop.int_to_binary(resistance_outcome[0])
        key2 = 'resistance ' + pop.int_to_binary(resistance_outcome[1])

        km_data_t[key1] = gen1_resistance_times
        km_data_t[key2] = gen2_resistance_times
        
        km_data[str(k_abs_t)] = km_data_t

    return km_data 