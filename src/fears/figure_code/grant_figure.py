from fears.classes.population_class import Population
import matplotlib.pyplot as plt
from fears.utils import plotter, results_manager
import numpy as np
#%%
np.random.seed(12)

init_counts = np.zeros(16)
init_counts[0] = 10**6

n_sims=1
options = {'n_sims':n_sims,
            'mut_rate':10**-5,
#            'plot':True,
            # 'death_rate':0.016,
            'k_abs':0.0008,
            'k_elim':0.0005,
            'max_dose':1000,
#            # 'death_rate':0,
            'constant_pop':False,
            'curve_type':'pharm',
            'timestep_scale':2,
            'n_timestep':4000,
            'plot':False,
            'carrying_cap':True,
            'max_cells':10**7,
            'death_rate':0.15,
            'init_counts':init_counts
            }

p_continuous = Population(**options)
p_discrete = Population(digital_seascape=True,mic_estimate = 0.2,**options)

p_continuous.simulate()
p_discrete.simulate()

#%%
fontsize=11
drug_kwargs = {'color':'black',
               'linewidth':2}

color_kwargs = {'palette':'tab20',
                'style':'solid'}
# color_kwargs={}

fig_digital, ax_digital = plt.subplots(2,1,figsize=(3.5,5))
t, ax_digital[0] = plotter.plot_fitness_curves(p_discrete,ax=ax_digital[0],
                                            show_legend=False,
                                            linewidth=2,
                                            labelsize=fontsize,
                                            color_kwargs=color_kwargs)

c = p_discrete.counts/np.max(p_discrete.counts)
p_discrete.drug_log_scale=False
dc = p_discrete.drug_curve/10**3

label = 'Drug Concentration \n($\u03BC$M x $10^{3}$)'
ax_digital[1],t = plotter.plot_timecourse_to_axes(p_discrete,
                                                c,
                                                ax_digital[1],
                                                drug_curve=dc,
                                                linewidth=2,
                                                drug_kwargs=drug_kwargs,
                                                labelsize=fontsize,
                                                drug_ax_sci_notation=False,
                                                drug_curve_label=label,
                                                color_kwargs=color_kwargs)

fig_sea, ax_sea = plt.subplots(2,1,figsize=(3.5,5))
t,ax_sea[0] = plotter.plot_fitness_curves(p_continuous,ax=ax_sea[0],
                                            show_legend=False,
                                            linewidth=2,
                                            labelsize=fontsize,
                                            color_kwargs=color_kwargs)

c = p_continuous.counts/np.max(p_continuous.counts)
p_continuous.drug_log_scale=False
dc = p_continuous.drug_curve/10**3

label = 'Drug Concentration \n($\u03BC$M x $10^{3}$)'
ax_sea[1],drug_ax = plotter.plot_timecourse_to_axes(p_continuous,
                                                c,
                                                ax_sea[1],
                                                drug_curve=dc,
                                                linewidth=2,
                                                drug_kwargs=drug_kwargs,
                                                labelsize=fontsize,
                                                drug_ax_sci_notation=False,
                                                drug_curve_label=label,
                                                color_kwargs=color_kwargs)

pos = ax_sea[1].get_position()
left = pos.x0
bottom = pos.y0 - 0.07
height = pos.height
width = pos.width
ax_sea[1].set_position([left, bottom, width, height])
ax_digital[1].set_position([left, bottom, width, height])
ax_sea[1].set_ylabel('Proportion of max cell count')
ax_digital[1].set_ylabel('Proportion of max cell count')
ax_sea[1].set_xlabel('Days')
ax_digital[1].set_xlabel('Days')
ax_sea[0].legend(loc=[1,0.05],ncol=2,frameon=False)
drug_ax.legend(loc=[1,1.3],frameon=False)
#%%

results_manager.save_fig(fig_digital,'discrete_seascape.png')
results_manager.save_fig(fig_sea,'continuous_seascape.png')