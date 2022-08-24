from fears.population import Population

p1 = Population(fitness_data='estimate')

p1.plot_fitness_curves()

p2 = Population(fitness_data='estimate',hc_estimate='joint')

p2.plot_fitness_curves()