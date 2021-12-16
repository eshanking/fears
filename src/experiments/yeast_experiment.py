# yeast experiment
from population_class import Population
import numpy as np

doubling_time = 1.5 # hours
death_rate = 0.016 # hours^-1
mut_rate = 10**-9
carrying_cap = True
max_cells = 10**11
n_timestep = 1000

init_counts = np.zeros(16)
init_counts[0] = 10**5

k_abs = 0.95 # hours^-1
k_elim = 0.00839 # hours^-1
max_dose = 1100 # uM??
curve_type = 'pulsed'
dose_schedule = 24
prob_drop = 0
pad_right = True

n_sims = 100
timestep_scale = 2


p1 = Population(doubling_time = doubling_time,
                death_rate = death_rate,
                mut_rate = mut_rate,
                k_abs = k_abs,
                k_elim = k_elim,
                max_dose = max_dose,
                curve_type = curve_type,
                n_sims = n_sims,
                dose_schedule = dose_schedule,
                prob_drop = prob_drop,
                carrying_cap = carrying_cap,
                max_cells = max_cells,
                pad_right = pad_right,
                init_counts=init_counts,
                n_timestep = n_timestep,
                constant_pop=False,
                timestep_scale=timestep_scale,
                plot=False)

c, n_survive = p1.simulate()
