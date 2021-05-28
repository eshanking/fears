from experiment_class_raw import Experiment
import numpy as np
import time
# import numpy as np

max_doses = [100]
curve_types = ['pharm']
experiment_type = 'rate-survival'
n_sims = 100
# slopes = [.0005,0.005,0.05]
# slopes = np.linspace(0.0005,0.0105,5)
# slopes = np.linspace(0.0005,0.0045,5)
slopes = np.array([0.0002,0.0004,0.0006,0.0008,0.001])
# slopes = np.linspace(0.0001,0.01,num=3)

# scale = 5 # scale time by two
mut_rate = 0.00005
death_rate = 0.3
# mut_rate = mut_rate/scale
# death_rate = death_rate/scale

options = {'mut_rate':mut_rate,
           'death_rate':death_rate,
           # 'x_lim':100,
           # 'y_lim':2*10**4,
           # 'counts_log_scale':True, 
            'plot':False,
            'k_elim':0.00009,
            'n_gen':2000,
            'ic50_path':'cycloguanil_ic50.csv'}

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