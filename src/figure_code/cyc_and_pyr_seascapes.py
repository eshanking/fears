from fears.utils import plotter, results_manager
from fears.classes.population_class import Population
import matplotlib.pyplot as plt

fig,ax = plt.subplots(nrows=2,ncols=1,sharey=False,figsize=(4,6))

p_pyr = Population(ic50_data='pyrimethamine_ic50.csv')
p_cyc = Population(ic50_data='cycloguanil_ic50.csv')

t,ax[0] = plotter.plot_fitness_curves(p_pyr,ax=ax[0],show_legend=False,
                                    linewidth=2)
t,ax[1] = plotter.plot_fitness_curves(p_cyc,ax=ax[1],show_legend=False,
                                    linewidth=2)

ax[0].set_ylabel('Growth rate ($hr^{-1}$)')
ax[1].set_ylabel('Growth rate ($hr^{-1}$)')

pad = -0.05

# ax[0].set_ylabel('% surviving')
pos1 = ax[1].get_position()
pos1.y0 = pos1.y0 + pad
pos1.y1 = pos1.y1 + pad
ax[1].set_position(pos1)

ax[1].legend(ncol=4,frameon=False,loc=(-0.1,-0.9),title='Genotypes')

results_manager.save_fig(fig,'pyr_cyc_seascapes.pdf')