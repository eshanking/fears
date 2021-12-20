from fears.classes.population_class import Population

"""
Note on mic estimate: the mic estimate parameter is a way to estimate the MIC
from the Hill equation. mic_estimate is the ratio of g/g_drugless (ratio of 
growth rate to max growth rate) at the MIC. For instance, if we assume the MIC
inhibits growth by 90%, then mic_estimate = 0.1.
"""

options = {'ic50_data':'pyrimethamine_ic50.csv'}
p_continuous = Population(**options) # population with continuous seascape
p_discrete = Population(digital_seascape=True,mic_estimate = 0.2,
                        **options) # population with digital/discrete seascape

p_continuous.plot_fitness_curves() # plot the fitness curves
p_discrete.plot_fitness_curves()

# generate fitness given a genotype and drug concentration
genotype = 5
conc = 100 # uM
f = p_continuous.gen_fitness(genotype,conc,p_continuous.drugless_rates,
                             p_continuous.ic50)

# generate the fitness landscape vector for a given drug concentration
fl = p_discrete.gen_fit_land(conc)