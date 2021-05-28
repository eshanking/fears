from population_class import Population
import numpy as np

# Generate a population, p, from fitness landscape data representing 32 alleles
p = Population(fitness_data='manual',landscape_path='test_fitness_data.txt',n_sims=10)

# counts = p.run_abm_v2() # this represents 1 simulation
# p.plot_timecourse(counts_t=counts)

# another option
# counts, n_survive = p.simulate() 
# this runs the simulation n_sims times (default 1) and averages the result
# it also stores counts as p.counts and outputs the result as well as the number of simulations that survived
# if p.plot is True, then p.simulate() will automatically plot the result of every simulation

# increase mutation rate and initial population size:
init_counts = np.zeros(32)    
init_counts[0] = 10**5

p = Population(fitness_data='manual',
               landscape_path='test_fitness_data.csv',
               n_sims=1,
               mut_rate=0.01,
               init_counts=init_counts,
               n_gen=2000)

counts, n_survive = p.simulate() 