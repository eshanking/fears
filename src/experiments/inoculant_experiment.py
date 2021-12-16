from experiment_class import Experiment

# experiment_type = 'inoculant-survival'
# options = {'n_gen':1000,'v2':True,'max_dose':200,'max_cells':10**7}
# inoculants = [10**3,10**4,10**5,10**6]
# curve_types=['constant','linear','pharm']
# n_sims = 10

# e = Experiment(experiment_type = experiment_type,
#                inoculants=inoculants,
#                curve_types=curve_types,
#                n_sims=n_sims,
#                population_options=options)
# e.run_experiment()
# e.plot_barchart()

# experiment_type = 'inoculant-survival'
# options = {'n_gen':1000,'v2':True,'max_dose':400,'max_cells':10**7}
# inoculants = [10**3,10**4,10**5,10**6]
# curve_types=['constant','linear','pharm']
# n_sims = 10

# e = Experiment(experiment_type = experiment_type,
#                inoculants=inoculants,
#                curve_types=curve_types,
#                n_sims=n_sims,
#                population_options=options)
# e.run_experiment()
# e.plot_barchart()

experiment_type = 'inoculant-survival'
options = {'n_gen':1000,'v2':True,'max_dose':800,'max_cells':10**7}
inoculants = [10**3,10**4,10**5,10**6]
curve_types=['constant','linear','pharm']
n_sims = 10

e = Experiment(experiment_type = experiment_type,
               inoculants=inoculants,
               curve_types=curve_types,
               n_sims=n_sims,
               population_options=options)
e.run_experiment()
e.plot_barchart()
