from fears.classes.experiment_class import Experiment
import numpy as np

np.random.seed(10937)

init_counts = np.zeros(16)
init_counts[0] = 10**10

transition_times = [800,1100]
options = {'doubling_time':1.5,
           'death_rate':0,
           'mut_rate':10**-6,
           'n_sims':10,
           # 'carrying_cap':False,
           'max_cells':10**9,
           'n_timestep':2000,
           'init_counts':init_counts,
           'timestep_scale':15,
           'plot':False,
           'ic50_data':'cycloguanil_ic50.csv',
           'drug_log_scale':True,
           'counts_log_scale':False}

e = Experiment(debug=False,
               experiment_type='ramp_up_down',
               ramp = 200,
               transition_times=transition_times,
               population_options=options)
e.run_experiment()