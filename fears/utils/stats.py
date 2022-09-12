import os
import numpy as np
from fears.utils import results_manager
import pickle
import lifelines

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
                extinction_time(pop,data,thresh=1)


            gen1_resistance_obs[k],gen1_resistance_times[k] = \
                resistance_time(pop,data,resistance_outcome[0],thresh=.1)
    
            gen2_resistance_obs[k],gen2_resistance_times[k] = \
                resistance_time(pop,data,resistance_outcome[1],thresh=.1)
                
            k+=1
            
        km_data_t['survival'] = death_event_times
        if type(resistance_outcome[0]) == list:
            key1 = 'resistance' + str(resistance_outcome[0])
        else:
            key1 = 'resistance ' + pop.int_to_binary(resistance_outcome[0])
        if type(resistance_outcome[1]) == list:
            key2 = 'resistance' + str(resistance_outcome[1])
        else:
            key2 = 'resistance ' + pop.int_to_binary(resistance_outcome[1])

        km_data_t[key1] = gen1_resistance_times
        km_data_t[key2] = gen2_resistance_times
        
        km_data[str(k_abs_t)] = km_data_t

    return km_data 

def gen_neighbors(pop,genotype):
    mut = range(pop.n_allele)
    neighbors = [genotype ^ (1 << m) for m in mut]

    return neighbors

def extinction_time(pop,counts,thresh=0):
    
    if len(counts.shape) > 1:
        c = np.sum(counts,axis=1)
    else:
        c = counts

    e = np.argwhere(c<=thresh)
    if len(e) == 0:
        event_obs = 0
        event_time = len(c)
    else:
        event_obs = 1
        event_time = e[0]
    
    timestep_scale = pop.timestep_scale
    event_time = event_time*timestep_scale
    
    return event_obs, event_time

def resistance_time(pop,counts,genotype,thresh=0.01):
    
    if thresh < 1:
        thresh = thresh*pop.carrying_cap
        
    if type(genotype) == list:
        
        times = []

        for g in genotype:
            if len(counts.shape) > 1:
                c = counts[:,g]
            else:
                c = counts
                
            e = np.argwhere(c>thresh)
            if len(e) == 0:
                # event_obs = 0
                times.append(len(c))
            else:
                # event_obs = 1
                times.append(e[0])
        
        if np.min(times) == len(c):
            event_time = len(c)
            event_obs = 0
        else:
            event_time = np.min(times)
            event_obs = 1
    
    else:
        if len(counts.shape) > 1:
            c = counts[:,genotype]
        else:
            c = counts
            
        e = np.argwhere(c>thresh)
        if len(e) == 0:
            event_obs = 0
            event_time = len(c)
        else:
            event_obs = 1
            event_time = e[0]
        
    timestep_scale = pop.timestep_scale
    event_time = event_time*timestep_scale
    
    return event_obs, event_time

def log_rank_test(self,durations_A, durations_B, 
                    event_observed_A=None, event_observed_B=None):
    
    results = lifelines.statistics.logrank_test(durations_A, durations_B, 
                                        event_observed_A=event_observed_A,
                                        event_observed_B=event_observed_B)
    
    return results