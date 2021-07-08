from fears.classes.experiment_class import Experiment
import numpy as np
import time
# import numpy as np
np.random.seed(2021)

max_doses = [10**4]
curve_types = ['pharm']
experiment_type = 'rate-survival'
n_sims = 100
# slopes = [.0005,0.005,0.05]
# slopes = np.linspace(0.0005,0.0105,5)
# slopes = np.linspace(0.0005,0.0045,5)
# slopes = np.array([0.0002,0.0004,0.0006,0.0008,0.001])
# slopes = np.array([.05])*10**-3
slopes = np.array([0.05,0.1,0.2])*10**-3
# slopes = np.linspace(0.0001,0.01,num=3)

init_counts = np.zeros(16)
init_counts[0] = 10**11

options = {'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-9,
           'carrying_cap':True,
           'max_cells':10**11,
           'n_timestep':2000,
           'init_counts':init_counts,
           # 'k_abs':0.95,
           # 'k_elim':0.00839,
           'k_elim':0,
           # 'max_dose':400,
           'dose_schedule':24,
           'pad_right':False,
           'timestep_scale':2,
           'plot':False,
           'ic50_data':'cycloguanil_ic50.csv'
           }

e = Experiment(max_doses=max_doses,
               slopes=slopes,
               curve_types = curve_types,
               experiment_type = experiment_type,
               n_sims=n_sims,
               population_options=options,
               debug=False)

t = time.time()
e.run_experiment()
# e.plot_barchart()
elapsed = time.time() - t
print(str(round(elapsed)))
