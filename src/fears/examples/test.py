#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 14:19:39 2021

@author: kinge2
"""

from fears.classes.population_class import Population
import numpy as np
# np.random.seed(1000)

p = Population(fitness_data='random',n_allele=2)
p.plot_fitness_curves()

p.set_null_seascape(10**0)
p.plot_fitness_curves()