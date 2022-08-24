from fears.population import Population

p1 = Population(fitness_data='estimate')

p1.plot_fitness_curves()

p2 = Population(fitness_data='estimate',hc_estimate='joint')

p2.plot_fitness_curves()

print('Drugless growth rates (per second):')
print("Method 1: " + str(p1.drugless_rates))
print("Method 2: " + str(p2.drugless_rates))

print('\n')
print('IC50s (log(ug/mL)):')
print("Method 1: " + str(p1.ic50))
print("Method 2: " + str(p2.ic50))

print('\n')
print('Hill coefficients:')
print("Method 1, genotype 0000: " + str(p1.seascape_lib['0']['hill_coeff']))
print("Method 2: " + str(p2.seascape_lib['0']['hill_coeff']))