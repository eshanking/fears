from fears.utils import results_manager
from fears.classes.experiment_class import Experiment
import numpy as np

e = Experiment()

p = e.populations[0]
p.death_rate = 0.2

fig = e.calculate_msw(0)

results_manager.save_fig(fig,'msw_0.pdf')