from fears.classes.population_class import Population

k_abs = 0.0001
k_elim = 0.000001
max_dose = 10**5
mut_rate = 10**-8
max_cells = 10**10
timestep_scale = 20
n_sims = 10
death_rate = 0

# ic50_data = 'pyrimethamine_ic50.csv'
ic50_data = 'cycloguanil_ic50.csv'

p_seascape = Population(ic50_data=ic50_data,
                        curve_type='pharm',
                        constant_pop=True,
                        k_abs=k_abs,
                        k_elim=k_elim,
                        max_dose=max_dose,
                        mut_rate=mut_rate,
                        max_cells=max_cells,
                        timestep_scale=timestep_scale,
                        n_sims=n_sims,
                        death_rate=death_rate)

p_landscape = Population(ic50_data=ic50_data,
                         curve_type='pharm',
                         static_landscape = True,
                         constant_pop=True,
                         k_abs=k_abs,
                         k_elim=k_elim,
                         max_dose=max_dose,
                         mut_rate=mut_rate,
                         max_cells=max_cells,
                         timestep_scale=timestep_scale,
                         n_sims=n_sims,
                         death_rate=death_rate)

p_seascape.simulate()
p_landscape.simulate()