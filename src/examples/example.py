from fears.classes.population_class import Population

# Non-zero death rate, non-constant population, with a carrying capacity
options = {
           'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-4,
           'max_cells':10**6,
           'plot':True,
           'k_abs':0.001,
           'k_elim':10**-4,
           'max_dose':10**1,
           'curve_type':'pharm',
           'timestep_scale':3,
           'n_sims':1,
           'debug':False
           }

p = Population(**options)
p.simulate()

# No death, constant population size, no carrying capacity
options = {
           'doubling_time':1.5,
           'death_rate':0,
           'constant_pop':True,
           'carrying_cap':False,
           'mut_rate':10**-4,
           'max_cells':10**6,
           'plot':True,
           'k_abs':0.001,
           'k_elim':10**-4,
           'max_dose':10**1,
           'curve_type':'pharm',
           'timestep_scale':3,
           'n_sims':1
           }

p = Population(**options)
p.simulate()

# Non-zero death, non-constant population, simulating a patient taking a daily
# drug dose regimen
options = {
           'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-4,
           'max_cells':10**6,
           'plot':True,
           'k_abs':0.1,
           'k_elim':10**-2,
           'max_dose':10**1,
           'curve_type':'pulsed',
           'dose_schedule':24,
           'timestep_scale':0.5,
           'n_sims':1,
           'debug':False
           }

p = Population(**options)
p.simulate()