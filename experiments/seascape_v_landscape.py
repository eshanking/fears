from fears.classes.experiment_class import Experiment
import numpy as np

def make_data():
    np.random.seed(109)
    
    init_counts = np.zeros(4)
    init_counts[0] = 10**10
    
    transition_times = [800,1600]
    options = {'doubling_time':1.5,
               'death_rate':0,
               'mut_rate':10**-9,
               'n_sims':10,
               # 'carrying_cap':False,
               'max_cells':10**11,
               'n_timestep':2400,
               'init_counts':init_counts,
               'timestep_scale':15,
               'plot':False,
               'plot_drug_curve':False,
               'drug_log_scale':True,
               'counts_log_scale':True,
               'max_dose':10**0}
    
    e = Experiment(debug=False,
                   experiment_type='ramp_up_down',
                   fitness_data = 'random',
                   n_allele=2,
                   ramp = 1,
                   transition_times=transition_times,
                   population_options=options,
                   # null_seascape_dose=10**-3,
                   null_seascape_dose=10**1,
                   first_dose=10**-3,
                   second_dose=10**1,
                   third_dose=10**5)
    
    e.run_experiment()
    
    return e

if __name__ == '__main__':
    make_data()