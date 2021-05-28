from experiment_class import Experiment

max_doses = [200]
curve_types = ['linear']
experiment_type = 'rate-survival'
n_sims = 10
slopes = [.1,.5,1.5]

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
            'plot':True}

e = Experiment(max_doses=max_doses,
               slopes=slopes,
               curve_types = curve_types,
               experiment_type = experiment_type,
               n_sims=n_sims,
               population_options=options)

e.run_experiment()
e.plot_barchart()
