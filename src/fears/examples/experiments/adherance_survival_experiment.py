from fears.classes.experiment_class import Experiment
import numpy as np

def make_data():
    np.random.seed(2021)
    
    init_counts = np.zeros(16)
    init_counts[0] = 10**5
    
    options = {'doubling_time':1.5,
               'death_rate':0.016,
               'mut_rate':10**-9,
               'carrying_cap':True,
               'max_cells':10**11,
               'n_timestep':1500,
               'init_counts':init_counts,
               'k_abs':0.95,
               'k_elim':0.00839,
               'max_dose':400,
               'dose_schedule':24,
               'pad_right':True,
               'timestep_scale':2,
               'plot':False,
               'ic50_data':'pyrimethamine_ic50.csv'}
    
    p = np.array([0.0,0.2,0.4,0.6,0.8])
    # p = np.array([0.6])
    # p = np.array([0.8])
    # n_sims = 500
    n_sims = 500
    
    experiment_type = 'drug-regimen'
    
    e = Experiment(experiment_type=experiment_type,
                   n_sims=n_sims,
                   prob_drops=p,
                   population_options = options,
                   debug=False)
    
    e.run_experiment()
    
    return e 

if __name__ == '__main__':
    make_data()