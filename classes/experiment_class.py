# Non-interactive class for running experiments and saving raw data.
# Does not produce images by default.

from fears.classes.population_class import Population
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
import warnings
import os
import time
import pickle
from fears.src import utils

class Experiment():
    
    # Initializer
    def __init__(self,
                 n_sims = 1,
                 curve_types = None,
                 max_doses = None,
                 low_dose=None,
                 mid_dose=None,
                 max_dose=None,
                 transition_times = None,
                 inoculants = None,
                 experiment_type = None,
                 prob_drops = None,
                 n_impulse=1,
                 population_options = {},
                 slopes=None,
                 debug = True): # debug = True -> no save
    
        self.root_path = str(utils.get_project_root())
        
        # list of allowed drug curve types
        allowed_types = ['linear',
                         'constant',
                         'heaviside',
                         'pharm',
                         'pulsed']
        
        # list of defined experiments
        allowed_experiments = ['inoculant-survival',
                               'dose-survival',
                               'drug-regimen',
                               'dose-entropy',
                               'rate-survival',
                               'bottleneck',
                               'ramp_up_down']
        
        if curve_types is not None:
            if not all(elem in allowed_types for elem in curve_types):
                raise Exception('One or more curve types is not recognized.\nAllowable types are: linear, constant, heaviside, pharm, pulsed.')
                
        if experiment_type is not None:
            if experiment_type not in allowed_experiments:
                raise Exception('Experiment type not recognized.\nAllowable types are inoculant-survival, dose-survival, drug-regimen, dose-entropy, and bottleneck.')
            
        # Curve type: linear, constant, heaviside, pharm, pulsed
        if curve_types is None:
            self.curve_types = ['constant']
        else:
            self.curve_types = curve_types
        
        if max_doses is None:
            self.max_doses = [1]
        else:
            self.max_doses = max_doses
            
        if inoculants is None:
            self.inoculants = [0]
        else:
            self.inoculants = inoculants
            
        self.n_sims = n_sims
            
        # Common options that are applied to each population
        self.population_options = population_options
        
        # initialize all populations
        self.populations = []
        
        # initialize a list of figures for saving
        self.figures = []
        
        if experiment_type is None:
            self.experiment_type = 'dose-survival'
            warnings.warn('No experiment type given - set to dose-survival by default.')
        else:
            self.experiment_type = experiment_type
        
        # if experiment_type == 'dose-survival' and len(inoculants) > 1:
        #     raise Exception('The experiment type is set to dose-survival (default), but more than one inoculant is given.')
        # elif experiment_type == 'inoculant-survival' and len(self.max_doses) > 1:
        #     # print('here')
        #     raise Exception('The experiment type is set to inoculant-survival, but more than one max dose is given.')
        
        # if experiment_type == 'inoculant-survival' and inoculants is None:
        #     raise Exception('The experiment type is set to inoculant-survival, but no inoculants are given.')
        
        if self.experiment_type == 'ramp_up_down':
            self.p_landscape = Population(consant_pop = True,
                                          carrying_cap = False,
                                          static_landscape = True,
                                          **self.population_options)
            self.p_seascape = Population(consant_pop = True,
                                          carrying_cap = False,
                                          **self.population_options)
            self.set_ramp_ud(self.p_landscape)
            self.set_ramp_ud(self.p_seascape)
        
        if self.experiment_type == 'bottleneck':
            for dc in self.duty_cycles:
                self.populations.append(Population(curve_type='on_off',
                                                   duty_cycle=dc,
                                                   n_sims = self.n_sims,
                                                   **self.population_options))
            # p_bottleneck = Population(max_dose = max_doses,
            #                              n_sims = self.n_sims,
            #                              curve_type = 'bottleneck'
            #                              **self.population_options)
            # self.populations.append(p_bottleneck)
            
            # p_no_bottleneck = Population(max_dose = max_doses,
            #                  n_sims = self.n_sims,
            #                  curve_type = 'constant'
            #                  **self.population_options)
            # self.populations.append(p_no_bottleneck)
        
        if self.experiment_type == 'dose-survival':
            for curve_type in self.curve_types:
                for max_dose in self.max_doses:
                    fig_title = 'Max dose = ' + str(max_dose) + ', curve type = ' + curve_type
                    self.populations.append(Population(curve_type=curve_type,
                                                      n_sims = self.n_sims,
                                                      max_dose = max_dose,
                                                      fig_title = fig_title,
                                                      **self.population_options))
                    
            self.n_survive = np.zeros([len(self.curve_types),len(self.max_doses)])
            self.perc_survive = np.zeros([len(self.curve_types),len(self.max_doses)])
                    
        elif self.experiment_type == 'inoculant-survival':
            for curve_type in self.curve_types:
                for inoculant in self.inoculants:
                    fig_title = 'Inoculant = ' + str(inoculant) + ', curve type = ' + curve_type
                    
                    init_counts = np.zeros(16)
                    init_counts[0] = inoculant
                    
                    self.populations.append(Population(curve_type=curve_type,
                                                      n_sims = self.n_sims,
                                                      fig_title = fig_title,
                                                      init_counts=init_counts,
                                                      **self.population_options))
                    
            self.n_survive = np.zeros([len(self.curve_types),len(self.inoculants)])
            self.perc_survive = np.zeros([len(self.curve_types),len(self.inoculants)])
            
        elif self.experiment_type == 'drug-regimen':
            
            self.prob_drops = prob_drops
            
            for prob_drop in self.prob_drops:
                curve_type = 'pulsed'
                self.populations.append(Population(curve_type=curve_type,
                                                   prob_drop=prob_drop,
                                                   # n_impulse = 1,
                                                   n_sims = 1,
                                                   # fig_title = fig_title,
                                                   # init_counts=init_counts,
                                            
                                                   **self.population_options))
            self.n_survive = np.zeros([len(self.populations)])
            
        elif self.experiment_type == 'dose-entropy':
            for dose in self.max_doses:
                self.populations.append(Population(max_dose = dose,
                                                   curve_type=self.curve_types[0]))
            self.entropy_results = pd.DataFrame(columns=[]) # will become a dataframe later
            
        elif self.experiment_type == 'rate-survival':
            # if the curve type is 'pharm' then slope will be interpreted as k_abs
            self.slopes = slopes
            for slope in self.slopes:
                if curve_types[0] == 'pharm':
                    self.populations.append(Population(max_dose=self.max_doses[0],
                                                        k_abs=slope,
                                                        curve_type='pharm',
                                                        n_sims=1,
                                                        **self.population_options))
                else:
                    self.populations.append(Population(max_dose=self.max_doses[0],
                                                        slope=slope,
                                                        curve_type='linear',
                                                        n_sims=1,
                                                        **self.population_options))                        
                    
            # self.rate_survival_results = pd.DataFrame(columns=[])
            
        # generate new save folder
                
        if not debug:
            num = 0
            num_str = str(num).zfill(4)
            
            date_str = time.strftime('%m%d%Y',time.localtime())
            
            save_folder = self.root_path + os.sep + 'results' + os.sep + 'results_' + date_str + '_' + num_str
            
            # save_folder = os.getcwd() + '//results_' + date_str + '_' + num_str
            
            while(os.path.exists(save_folder)):
                num += 1
                num_str = str(num).zfill(4)
                save_folder = self.root_path + os.sep + 'results' + os.sep + 'results_' + date_str + '_' + num_str
            os.mkdir(save_folder) 
            
            self.experiment_info_path = self.root_path + os.sep + 'results' + os.sep + 'experiment_info_' + date_str + '_' + num_str + '.p'
            
            # self.experiment_info_path = experiment_info_path
            self.results_path = save_folder
        
            # save the experiment information
            
            pickle.dump(self, open(self.experiment_info_path,"wb"))
        
        # self.n_survive = np.zeros([len(self.curve_types),len(self.max_doses)])
        # self.perc_survive = np.zeros([len(self.curve_types),len(self.max_doses)])
###############################################################################
    # Methods for running experiments
    
    # run experiment and save results
    def run_experiment(self):
            
        n_doses = len(self.max_doses)
        n_curves = len(self.curve_types)
        n_inoc = len(self.inoculants)
        
        # pbar = tqdm(total = n_curves*n_doses) # progress bar
        
        # Loop through each population, execute simulations, and store survival statistics
        
        if self.experiment_type == 'dose-survival':
            # pbar = tqdm(total = n_curves*n_doses) # progress bar
            for curve_number in range(n_curves):
                for dose_number in range(n_doses):
                    
                    exp_num = curve_number*n_doses + dose_number
                    pop = self.populations[exp_num] # extract population in list of population
                    c,n_survive_t = pop.simulate()
                    pop.plot_timecourse()
                    self.n_survive[curve_number,dose_number] = n_survive_t
                    # pbar.update()
            self.perc_survive = 100*self.n_survive/self.n_sims   
                 
        elif self.experiment_type == 'inoculant-survival':
            # pbar = tqdm(total = n_curves*n_inoc) # progress bar
            for curve_number in range(n_curves):
                for inoc_num in range(n_inoc):
                    
                    exp_num = curve_number*n_inoc + inoc_num
                    pop = self.populations[exp_num] # extract population in list of population
                    c,n_survive_t = pop.simulate()
                    pop.plot_timecourse()
                    self.n_survive[curve_number,inoc_num] = n_survive_t
                    # pbar.update()           
            self.perc_survive = 100*self.n_survive/self.n_sims
            
        elif self.experiment_type == 'drug-regimen':
            # pbar = tqdm(total=len(self.populations))
            # kk=0
            for p in self.populations:
                for i in range(self.n_sims):
                    # initialize new drug curve
                    p.drug_curve,u = p.gen_curves()
                    counts,n_survive = p.simulate()
                    drug = p.drug_curve
                    drug = np.array([drug])
                    drug = np.transpose(drug)
                    
                    regimen = self.compute_regimen(p, u)
                    regimen = np.array([regimen])
                    regimen = np.transpose(regimen)
                    regimen_t = np.zeros(drug.shape)
                    regimen_t[0:len(regimen)] = regimen
                    
                    counts = np.concatenate((counts,drug,regimen_t),axis=1)
  
                    save_folder = 'p_drop=' + str(p.prob_drop)
                    save_folder = save_folder.replace('.',',')
                    self.save_counts(counts,i,save_folder)
                # kk+=1
                # pbar.update()
                self.perc_survive = 100*self.n_survive/self.n_sims
            
        elif self.experiment_type == 'dose-entropy':
            # pbar = tqdm(total=len(self.populations)*self.n_sims)
            e_survived = []
            e_died = []
            for p in self.populations:
                
                for i in range(self.n_sims):
                    c,n_survive = p.simulate()
                    # e = max(p.entropy()) # compute max entropy
                    e_t = p.entropy()
                    e = max(e_t)
                    # e=1
                    # p.plot_timecourse()
                    
                    if n_survive == 1:
                        survive = 'survived' # survived
                        e_survived.append(e_t)
                    else:
                        survive = 'extinct' # died
                        e_died.append(e_t)      
                        
                    d = {'dose':[p.max_dose],
                         'survive condition':[survive],
                         'max entropy':[e]}
                    
                    entropy_results_t = pd.DataFrame(d)
                    self.entropy_results = self.entropy_results.append(entropy_results_t)
                    # pbar.update()
        
        elif self.experiment_type == 'rate-survival':
            # pbar = tqdm(total=len(self.populations))
            
            for p in self.populations:
                for n in range(self.n_sims):
                    counts,n_survive = p.simulate()
                    
                    drug = p.drug_curve
                    drug = np.array([drug])
                    drug = np.transpose(drug)
                    counts = np.concatenate((counts,drug),axis=1)
                    
                    if self.curve_types[0] == 'pharm':
                        save_folder = 'k_abs=' + str(p.k_abs)
                        save_folder.replace('.','pnt')
                    else:
                        save_folder = 'slope=' + str(p.slope)
                        save_folder.replace('.','pnt')
                    self.save_counts(counts,n,save_folder)
        
        
                    
                # fig_savename = 'slope = ' + str(p.slope)
                # self.figures = self.figures.append(fig)
                # pbar.update()
            # self.rate_survival_results.index = np.arange(len(self.rate_survival_results))
                
        # pbar.close() # close progress bar
  
    # save counts as a csv in the given subfolder with the label 'num'
    def save_counts(self,counts,num,save_folder,prefix='sim_'):
        
        # check if the desired save folder exists. If not, create it
        folder_path = self.results_path + os.sep + save_folder
        if os.path.exists(folder_path) != True:
            os.mkdir(folder_path)
            
        num = str(num).zfill(4)
        savename = self.results_path + os.sep + save_folder + os.sep + prefix + num + '.csv'
        np.savetxt(savename, counts, delimiter=",")
        return
    
    def compute_regimen(self,p,u):
        gap = int(p.dose_schedule/p.timestep_scale)
        n_impulse = int(np.ceil(p.n_timestep/gap))
        regimen = np.zeros(n_impulse)
        
        for i in range(n_impulse):
            if u[(i)*gap] == 1:
                regimen[i] = 1
                
        return regimen
    
    def set_ramp_ud(self,p):
        
        drug_curve = np.zeros()
        
        return drug_curve
        
###############################################################################
