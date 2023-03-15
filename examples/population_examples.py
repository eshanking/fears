from fears.population import Population

#%% Default data from Ogbunugafor 2016 (engineered S. cerevisiae subject to pyrimethamine)
p = Population(drug_units = 'uM')

p.ic50 = [i + 6 for i in p.ic50]

p.plot_fitness_curves()
# p.print_params()
#%% Default data from King 2022 (engineered E. coli subject to cefotaxime)

p = Population(fitness_data='from_file') # no file path -> default data
p.plot_fitness_curves()
p.print_params()