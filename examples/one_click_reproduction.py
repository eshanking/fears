# one click code reproduction
from fears.figure_code import pyr_seascape
from fears.experiments import rate_survival_experiment_pharm
from fears.experiments import adherance_survival_experiment
from fears.experiments import seascape_v_landscape

#%% Figure 1

pyr_seascape.make_figure()

#%% Experiment data

# rate of change vs survival
rate_survival_experiment_pharm.make_data()

# adherance vs survival
adherance_survival_experiment.make_data()

# seascape vs landscape
seascape_v_landscape.make_data()

#%% Make figures