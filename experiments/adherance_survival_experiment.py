from fears.classes.experiment_class_raw import Experiment
import numpy as np

init_counts = np.zeros(16)
init_counts[0] = 10**5

options = {'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-9,
           'carrying_cap':True,
           'max_cells':10**11,
           'n_timestep':1000,
           'init_counts':init_counts,
           'k_abs':0.95,
           'k_elim':0.00839,
           'max_dose':1100,
           'dose_schedule':24,
           'pad_right':True,
           'timestep_scale':2,
           'plot':False}

p = np.array([0.2])
n_sims = 1000

experiment_type = 'drug-regimen'

e = Experiment(experiment_type=experiment_type,
               n_sims=n_sims,
               prob_drops=p,
               population_options = options,
               debug=False)

e.run_experiment()
