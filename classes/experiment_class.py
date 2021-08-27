from fears.classes.population_class import Population
from fears.utils import plotter, dir_manager
import numpy as np
import pandas as pd
# import warnings
import os
import time
import pickle
import lifelines

class Experiment():

    # Initializer
    def __init__(self,
                 n_sims = 1,
                 curve_types = None,
                 max_doses = None,
                 first_dose = None,
                 third_dose = None,
                 second_dose = None,
                 ramp = 100,
                 fitness_data = 'generate',
                 n_allele = 4,
                 null_seascape_dose=None,
                 transition_times = None,
                 inoculants = None,
                 experiment_type = None,
                 prob_drops = None,
                 n_impulse=1,
                 population_options = {},
                 slopes=None,
                 passage=False,
                 passage_time = 48,
                 debug = True): # debug = True -> no save
    
        self.root_path = str(dir_manager.get_project_root())
        
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
            # warnings.warn('No experiment type given - set to dose-survival by default.')
        else:
            self.experiment_type = experiment_type
        
        # if experiment_type == 'dose-survival' and len(inoculants) > 1:
        #     raise Exception('The experiment type is set to dose-survival (default), but more than one inoculant is given.')
        # elif experiment_type == 'inoculant-survival' and len(self.max_doses) > 1:
        #     # print('here')
        #     raise Exception('The experiment type is set to inoculant-survival, but more than one max dose is given.')
        
        # if experiment_type == 'inoculant-survival' and inoculants is None:
        #     raise Exception('The experiment type is set to inoculant-survival, but no inoculants are given.')
        self.ramp = ramp
        if self.experiment_type == 'ramp_up_down':
            self.p_landscape = Population(constant_pop = True,
                                          carrying_cap = False,
                                          fitness_data=fitness_data,
                                          n_allele=n_allele,
                                          **self.population_options)
            self.p_seascape = Population(constant_pop = True,
                                         carrying_cap = False,
                                         fitness_data=fitness_data,
                                         n_allele=n_allele,
                                         **self.population_options)
            
            if fitness_data=='random':
                self.p_seascape.ic50 = self.p_landscape.ic50
                self.p_seascape.drugless_rates = self.p_landscape.drugless_rates
                self.p_landscape.set_null_seascape(null_seascape_dose)
            self.null_seascape_dose = null_seascape_dose
            
            if second_dose is None:
                self.second_dose = 10**5
            else:
                self.second_dose=second_dose
            if third_dose is None:
                self.third_dose = 10**1
            else:
                self.third_dose=third_dose
            if first_dose is None:
                self.first_dose = 10**-2
            else:
                self.first_dose = first_dose
            if transition_times is None:
                n_timestep = self.p_landscape.n_timestep
                self.transition_times = [int(n_timestep/4),int(3*n_timestep/4)]
            else:
                self.transition_times = transition_times
            
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
                                                        passage=passage,
                                                        passage_time=passage_time,
                                                        **self.population_options))
                elif curve_types[0] == 'pulsed':
                    self.populations.append(Population(max_dose=self.max_doses[0],
                                                        k_abs=slope,
                                                        curve_type='pulsed',
                                                        n_sims=1,
                                                        passage=passage,
                                                        passage_time=passage_time,
                                                        **self.population_options))
                else:
                    self.populations.append(Population(max_dose=self.max_doses[0],
                                                        slope=slope,
                                                        curve_type='linear',
                                                        n_sims=1,
                                                        passage=passage,
                                                        passage_time=passage_time,
                                                        **self.population_options))                        
                    
            # self.rate_survival_results = pd.DataFrame(columns=[])
            
        # generate new save folder
        
        self.debug=debug
        
        if not debug:
            
            num = 0
            num_str = str(num).zfill(4)
            
            date_str = time.strftime('%m%d%Y',time.localtime())
            
            save_folder = self.root_path + os.sep + 'results' + os.sep + 'results_' + date_str + '_' + num_str
            
            # save_folder = os.getcwd() + '//results_' + date_str + '_' + num_str
            if not (os.path.exists(self.root_path + os.sep + 'results')):
                os.mkdir(self.root_path + os.sep + 'results')
                
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
        
        elif self.experiment_type == 'ramp_up_down':
            counts_landscape, ft = self.p_landscape.simulate()
            drug_curve = self.p_landscape.drug_curve
            drug_curve= np.array([drug_curve])
            drug_curve = np.transpose(drug_curve)
            counts_seascape, ft = self.p_seascape.simulate()
            
            if not self.debug:
                savedata = np.concatenate((counts_landscape,drug_curve),axis=1)
                self.save_counts(savedata, num=None, save_folder=None,prefix = 'landscape_counts')
                
                savedata = np.concatenate((counts_seascape,drug_curve),axis=1)
                self.save_counts(savedata, num=None, save_folder=None,prefix = 'seascape_counts')
        
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
                    if not self.debug:
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
                    
                    if self.debug is False:
                        if (self.curve_types[0] == 'pharm' or 
                            self.curve_types[0] == 'pulsed'):
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
        if save_folder is None:
            save_folder = ''
        folder_path = self.results_path + os.sep + save_folder
        if os.path.exists(folder_path) != True:
            os.mkdir(folder_path)
        
        if num is None:
            num = ''
        else:
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
        
        n_timestep = p.n_timestep
        drug_curve = np.zeros(n_timestep)
        slope = (self.second_dose-self.first_dose)/self.ramp
        buffer = self.ramp/2
        
        times = [self.transition_times[0]-buffer,
                 self.transition_times[0]+buffer,
                 self.transition_times[1]-buffer,
                 self.transition_times[1]+buffer]
        
        for t in range(n_timestep):
            if t<times[0]:
                drug_curve[t] = self.first_dose
            elif t<times[1]:
                drug_curve[t] = drug_curve[t-1]+slope
            elif t<times[2]:
                drug_curve[t] = self.second_dose
            elif (t<times[3] and 
                  drug_curve[t-1]+slope<self.third_dose):
                drug_curve[t] = drug_curve[t-1]+slope
            else:
                drug_curve[t] = self.third_dose
                
        p.drug_curve = drug_curve
        
        
        return drug_curve
    
    def calculate_msw(self,genotype,pop=None,save=False):
        if pop is None:
            pop=self.populations[0]
        
        # calculate neighbors in bit-string network model
        neighbors = self.gen_neighbors(pop,genotype)
        genotypes = [genotype] + neighbors
        
        powers = np.linspace(-3,5,40)
        conc = np.power(10*np.ones(powers.shape[0]),powers)
        
        fitness_curves = np.zeros((len(powers),len(genotypes)))
        
        drugless_rates = pop.drugless_rates
        ic50 = pop.ic50
        
        for j in range(len(genotypes)):
            for i in range(len(conc)):
                
                f = pop.gen_fitness(genotypes[j],
                                    conc[i],
                                    drugless_rates,
                                    ic50)
                r_0 = f-pop.death_rate
                fitness_curves[i,j] = r_0
        fig = plotter.plot_msw(pop,fitness_curves,conc,genotypes,save=save)
        return fig
        
###############################################################################
# Helper methods
    def gen_neighbors(self,pop,genotype):
        mut = range(pop.n_allele)
        neighbors = [genotype ^ (1 << m) for m in mut]

        return neighbors
    
    def extinction_time(self,pop,counts,thresh=1):
        
        if len(counts.shape) > 1:
            c = np.sum(counts,axis=1)
        else:
            c = counts
        e = np.argwhere(c<thresh)
        if len(e) == 0:
            event_obs = 0
            event_time = len(c)
        else:
            event_obs = 1
            event_time = e[0]
        
        timestep_scale = pop.timestep_scale
        event_time = event_time*timestep_scale
        
        return event_obs, event_time
    
    def resistance_time(self,pop,counts,genotype,thresh=0.01):
        
        if len(counts.shape) > 1:
            c = counts[:,genotype]
        else:
            c = counts
        
        if thresh < 1:
            thresh = thresh*pop.max_cells
            
        e = np.argwhere(c>thresh)
        if len(e) == 0:
            event_obs = 0
            event_time = len(c)
        else:
            event_obs = 1
            event_time = e[0]
        
        timestep_scale = pop.timestep_scale
        event_time = event_time*timestep_scale
        
        return event_obs, event_time
    
    def log_rank_test(self,durations_A, durations_B, 
                      event_observed_A=None, event_observed_B=None):
        
        results = lifelines.statistics.logrank_test(durations_A, durations_B, 
                                          event_observed_A=event_observed_A,
                                          event_observed_B=event_observed_B)
        
        return results
    