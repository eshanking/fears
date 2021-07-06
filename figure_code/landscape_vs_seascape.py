from fears.classes.population_class import Population

p = Population(static_topology=True,ic50_data='cycloguanil_ic50.csv',static_topo_dose=10**0)
p.plot_fitness_curves(fig_title = 'Landscape',save=True,savename='static_topology.pdf')

p = Population(static_topology=False,ic50_data='cycloguanil_ic50.csv')
p.plot_fitness_curves(fig_title='Seascape',save=True,savename='seascape.pdf')