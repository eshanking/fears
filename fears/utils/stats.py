import os
import numpy as np
from fears.utils import results_manager
import pickle
import lifelines

def survival_proportion(pop,data):
    """Computes survival fraction of populations in an experiment for lhs analysis 

    Args:
        pop (population): Population class object
        data (list): list of population size counts over time

    Returns:
        float: proportion of simulations that survived
    """
    n_survived = 0
    for c in data:
        obs,time = extinction_time(pop,c,thresh=1)
        if obs == 0:
            n_survived += 1
    
    p_survived = n_survived/len(data)

    return p_survived

def generate_binary_strings(N, m):
    def generate_strings_helper(N, m, prefix):
        if N == 0:
            if m == 0:
                binary_str = ''.join(prefix)
                decimal_value = int(binary_str, 2)
                results.append(decimal_value)
            return
        if m > 0:
            generate_strings_helper(N - 1, m - 1, prefix + ['1'])
        generate_strings_helper(N - 1, m, prefix + ['0'])

    results = []
    generate_strings_helper(N, m, [])
    return results

def most_freq_genotype(exp=None,exp_info_path=None,mode='mode'):
    """Computes the proportion of the population that has a given number of mutations
    
    Returns a dictionary of dictionaries of K-M curves from the 
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

    max_idx_dict = {}
    
    pop = exp_info.populations[0]
    
    for exp in exp_folders:
    
        k_abs_t = exp[exp.find('=')+1:]
        k_abs_t = k_abs_t.replace(',','.')
        k_abs_t = float(k_abs_t)
        
        k_abs_t = round(k_abs_t,10)
        # print(f"{k_abs_t:.2e}")
        
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        
        if mode == 'mode':
            k=0

            data = np.zeros((pop.n_timestep,pop.n_genotype))

            while k < len(sim_files):

                sim = sim_files[k]
                sim = exp + os.sep + sim
                data_dict = results_manager.get_data(sim)

                data += data_dict['counts']
                k+=1
                
            max_idx = np.argmax(data,axis=1)
            max_idx = max_idx[data.sum(axis=1) > 0]
            max_idx_dict[str(k_abs_t)] = max_idx
            
        else:
            k=0

            data = np.zeros((pop.n_timestep,len(sim_files)))
            pop_size = np.zeros((pop.n_timestep))
            while k < len(sim_files):
                    
                sim = sim_files[k]
                sim = exp + os.sep + sim
                data_dict = results_manager.get_data(sim)

                pop_size += np.sum(data_dict['counts'],axis=1)

                data[:,k] = np.argmax(data_dict['counts'],axis=1)
                k+=1
            
            data = data[pop_size > 0,:]
            max_idx_dict[str(k_abs_t)] = data
        
    return max_idx_dict

def n_mut_curve(exp=None,exp_info_path=None,nmut=1):
    """Computes the proportion of the population that has a given number of mutations
    
    Returns a dictionary of dictionaries of K-M curves from the 
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

    genotypes = generate_binary_strings(pop.n_allele,nmut)
    
    prop_data = {}

    # c[:,0] = c[:,0] + np.ones(len(c[:,0]))

    # pop_size = np.sum(c,axis=1)

    # for i in range(c.shape[1]):
    #     c_t = c[:,i]
    #     c[:,i] = np.divide(c_t,pop_size)
    
    for exp in exp_folders:
    
        k_abs_t = exp[exp.find('=')+1:]
        k_abs_t = k_abs_t.replace(',','.')
        k_abs_t = float(k_abs_t)
        
        k_abs_t = round(k_abs_t,10)
        # print(f"{k_abs_t:.2e}")
        
        sim_files = os.listdir(path=exp)
        sim_files = sorted(sim_files)
        
        k=0

        data = np.zeros((pop.n_timestep,pop.n_genotype))

        while k < len(sim_files):

            sim = sim_files[k]
            sim = exp + os.sep + sim
            data_dict = results_manager.get_data(sim)

            data += data_dict['counts']
            k+=1
            
        pop_size = np.sum(data,axis=1)
        data = data[pop_size > 0]
        pop_size = pop_size[pop_size > 0]
        data = np.divide(data,pop_size[:,None])
        data = data[:,genotypes]
        data = np.sum(data,axis=1)

        prop_data[str(k_abs_t)] = data
        
    return prop_data

def km_curve(exp=None,exp_info_path=None,resistance_outcome=[14,15]):
    """Computes Kaplan-Meier data for km curve plotting
    
    Returns a dictionary of dictionaries of K-M curves from the 
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