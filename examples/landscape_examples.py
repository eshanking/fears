from fears.population import Population
from fears.utils import fitness,plotter

p = Population(fitness_data='from_file')

# plot the fitness seascape
fig,ax = p.plot_fitness_curves()

# plot some fitness landscapes at different drug concentrations
plotter.plot_landscape(p,10**-3)

plotter.plot_landscape(p,10**2)

# enumerate neighbors

print(p.gen_neighbors(0)) # genotype 0 (0000) neighbors

print(p.gen_neighbors(15)) # genotype 15 (1111) neighbors 

gen_15_neighbors = p.gen_neighbors(15)
gen_15_neighbors = [p.int_to_binary(g) for g in gen_15_neighbors]

print(gen_15_neighbors) # now in binary format