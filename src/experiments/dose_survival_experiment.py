from experiment_class import Experiment

max_doses = [.1,1,5]
curve_types = ['constant','linear','pharm']
experiment_type = 'dose-survival'
n_sims = 10

options = {'mut_rate':0.001,
           'death_rate':0.50}

e = Experiment(max_doses=max_doses,
               curve_types = curve_types,
               experiment_type = experiment_type,
               n_sims=n_sims,
               population_options=options)

e.run_experiment()
e.plot_barchart()
