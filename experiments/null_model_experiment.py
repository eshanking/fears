from fears.classes.population_class import Population
import numpy as np
import matplotlib.pyplot as plt

# k_abs = 0.0001
# k_elim = 0.000001
# max_dose = 10
# mut_rate = 10**-2
# max_cells = 10**6
# timestep_scale = 2
# n_timestep = 1000
# n_sims = 10
# death_rate = 0

# # ic50_data = 'pyrimethamine_ic50.csv'
# ic50_data = 'cycloguanil_ic50.csv'

# p_seascape = Population(ic50_data=ic50_data,
#                         curve_type='constant',
#                         constant_pop=True,
#                         carrying_cap=False,
#                         k_abs=k_abs,
#                         k_elim=k_elim,
#                         max_dose=max_dose,
#                         mut_rate=mut_rate,
#                         max_cells=max_cells,
#                         timestep_scale=timestep_scale,
#                         n_sims=n_sims,
#                         death_rate=death_rate,
#                         stop_condition=False,
#                         n_timestep=n_timestep)

# p_landscape = Population(ic50_data=ic50_data,
#                          curve_type='pharm',
#                          static_landscape = True,
#                          constant_pop=True,
#                          carrying_cap=False,
#                          k_abs=k_abs,
#                          k_elim=k_elim,
#                          max_dose=max_dose,
#                          mut_rate=mut_rate,
#                          max_cells=max_cells,
#                          timestep_scale=timestep_scale,
#                          n_sims=n_sims,
#                          death_rate=death_rate,
#                          stop_condition=True,
#                         n_timestep=n_timestep)

# p_seascape.simulate()
# p_landscape.simulate()

init_counts = np.zeros(16)
init_counts[0] = 10**5

options = {
           'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-9,
           'carrying_cap':False,
           'max_cells':10**11,
           'n_timestep':1000,
           'init_counts':init_counts,
           'k_abs':0.001,
           'k_elim':10**-5,
           'max_dose':10**3,
           'dose_schedule':24,
           'pad_right':True,
           'timestep_scale':2,
           'plot':True,
           'ic50_data':'cycloguanil_ic50.csv',
           'constant_pop':True,
           'stop_condition':False,
           'curve_type':'pharm',
           'n_sims':10,
           'plot':True
           }

p_seascape = Population(**options)

seascape_fix_times = np.array(p_seascape.simulate())*(p_seascape.timestep_scale/24)

p_landscape = Population(static_landscape=True,
                          **options)

landscape_fix_times = np.array(p_landscape.simulate())*(p_landscape.timestep_scale/24)

fig,ax = plt.subplots(figsize=(3,3))
ax.boxplot([seascape_fix_times,landscape_fix_times])

x = np.ones(len(seascape_fix_times))
ax.scatter(x,seascape_fix_times)

x = np.ones(len(landscape_fix_times))*2
ax.scatter(x,landscape_fix_times)

ax.set_ylabel('Fixation time (days)',fontsize=12)
ax.set_xticklabels(['Seascape','Landscape'],fontsize=12)

ax.set_ylim([10,80])

# plt.savefig('fixation_time.pdf',bbox_inches="tight")