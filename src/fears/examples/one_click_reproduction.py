# one click code reproduction
# Note: this will take several hours to complete

from fears.figure_code import pyr_seascape, example_timecourses, rate_of_change_km_fig, nonadherance_km_figure
from fears.experiments import rate_survival_experiment_pharm
from fears.experiments import adherance_survival_experiment
from fears.experiments import seascape_v_landscape

#%% Figure 1

pyr_seascape.make_figure()

#%% Experiment data

# rate of change vs survival
roc_exp = rate_survival_experiment_pharm.make_data()

# adherance vs survival
adh_exp = adherance_survival_experiment.make_data()

# seascape vs landscape
svl_exp = seascape_v_landscape.make_data()

#%% Make figures

#%% Figure 1 - pyrimethamine seascape

pyr_seascape.make_figure()

#%% Figure 2 - representative data

example_timecourses.make_figure(roc_exp,adh_exp)

#%% Figure 3 - rate of change K-M curves
# Also generates p-value tables

rate_of_change_km_fig.make_fig(roc_exp)

#%% Figure 4 - nonadherance K-M curves
# Also generates p-value tables

nonadherance_km_figure.make_fig(adh_exp)