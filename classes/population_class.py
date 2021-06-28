import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# from cycler import cycler
# import seaborn as sns
# import scipy as sp
import math
from fears.utils import plotter, pharm, fitness, dir_manager

class Population:
    
    """Population class: the fundamental object of FEArS.
    
    Contains basic information about the organism being simulated, the 
    environment, and the evolutionary algorithm.
    ...
    
    Attributes
    __________
    passage : bool
        if true, simulates passaging cells by reducing the population to 10% 
        at intervals given by drug_regimen.
    carrying_cap : bool
        if true, sets the carrying capacity to max_cells.
    curve_type : str
        determines the drug concentration curve. 
        Allowed types:
            constant: constant at max dose.
            linear: linear ramp to max_dose at a rate given by slope.
            heaviside: jumps drug concentration from min_dose to max_dose at a
                timestep given by h_step
            pharm: drug curve follows 1-compartment pharmacokinetic model. 
                k_abs: absorption coefficient
                k_elim: eliminiation coefficient
            pulsed: simulates a patient taking a drug dose at intervals given 
                by dose_schedule
            on_off: switches between 'on' (max_dose) and 'off' (min_dose) at
                intervals given by dose_schedule. On/off ratio set by 
                duty_cycle
    counts_log_scale : bool
        if true, plots the results on a log scale
    constant_pop : enforces a constant population size (max_cells)
    drugless_data : str
        filename for the drugless growth rates (ogbunugafor_ic50.csv by 
        default). Searches for files in the fears/data folder.
    death_rate : float
        death rate
    doubling_time : float
        average doubling time of the model organism (hours)
    dose_schedule : int
        timesteps between simulated doses
    drug_log_scale : bool
        if true, plots the drug concentration curve on a log scale
    drug_curve : array
        optional custom drug concentration curve. Overrides all other drug 
        concentration options.
    debug : bool
        if true, prints every 10th timestep
    duty_cycle : float
        on/off ratio of the on/off regimen.
    entropy_lim : list
        y axis limits for plotting the entropy curve
    fig_title : str
        optional figure title
    fitness_data : str
        generate: generate the fitness data based on the drugless growth rates
        and ic50.
        manual: get fitness data from self.landscape_data
    h_step: 
    
    """
###############################################################################    
    # Initializer
    def __init__(self,
                 bottleneck = False,
                 carrying_cap = True,
                 curve_type='constant', # drug concentration curve
                 counts_log_scale = False, # plot counts on log scale
                 constant_pop = False, # normalize to a constant population size
                 drugless_data = None, # file path for the drugless growth rates
                 death_rate = 0.15, 
                 doubling_time = 1, # average doubling time of model organism
                 dose_schedule=12, # dose every x hours
                 drug_log_scale = False, # plot the drug concentration curve on a log scale
                 drug_curve = None, # input a custom drug concentration curve
                 debug=False, # print the current time step
                 duty_cycle = None,
                 entropy_lim = None, # entropy plotting limits
                 fig_title = '',
                 fitness_data = 'generate', # 'generate' = generate fitness data using drugless growth rates, ic50, and drug concentration. 'manual' = input fitness landscape from csv
                 h_step = 500,
                 ic50_data = None, 
                 init_counts = None, # default is 10,000 wild type cells
                 k_elim = 0.001, # for modeling pharmacokinetics
                 k_abs = 0.07,
                 landscape_path = None, # path for custom fitness landscape
                 min_dose = 0,
                 mut_rate = 0.01, # mutation rate
                 max_cells = 10**6, # carrying capacity
                 max_dose = 1, 
                 n_timestep=1000, # number of generations
                 n_sims = 1, # number of simulations to average together
                 pad_right = False,
                 plot=True, # plot the result of simulate()
                 plot_entropy = False, # plot the entropy of the population over time underneath the timecourse
                 prob_drop=0, 
                 slope = None, 
                 static_topology = False,
                 static_topo_dose = 10**5,
                 stop_condition = False,
                 timestep_scale = 1,
                 x_lim = None, # plotting
                 y_lim = None # plotting
                 ):
                
        # Evolutionary parameters
        
        # Number of generations (time steps)
        
        if carrying_cap is False and constant_pop is False:
            print('\nWarning: no limit set on population size! Consider'
                          + ' enforcing a carrying capacity or constant'
                          + ' population.')
        
        self.n_timestep = n_timestep
        self.stop_condition = stop_condition
        self.max_cells = max_cells
        
        # model parameters
        self.mut_rate = mut_rate
        self.death_rate = death_rate
        self.doubling_time = doubling_time
        self.timestep_scale = timestep_scale # timestep_scale = 2 -> timestep = 2 hrs, etc
        
        self.carrying_cap = carrying_cap
        self.n_sims = n_sims # number of simulations to average together in self.simulate
        self.constant_pop = constant_pop
        self.debug = debug
        
        self.fitness_data = fitness_data
        
        self.bottleneck = bottleneck

        self.counts = np.zeros([self.n_timestep,16])
        self.counts_extinct = np.zeros([self.n_timestep,16])
        self.counts_survive = np.zeros([self.n_timestep,16])
        
        # Generate fitness data from IC50 and drugless growth rate data
        if fitness_data == 'generate':
            # Data paths
            if drugless_data is None:
                self.drugless_data = dir_manager.make_datapath_absolute('ogbunugafor_drugless.csv')
            else:
                self.drugless_data = dir_manager.make_datapath_absolute(drugless_data)
                
            if ic50_data is None:
                # self.ic50_data = "C:\\Users\\Eshan\\Documents\\python scripts\\theory division\\abm_variable_fitness\\data\\pyrimethamine_ic50.csv"
                # self.ic50_data = self.make_datapath_absolute('pyrimethamine_ic50.csv')
                self.ic50_data = dir_manager.make_datapath_absolute('pyrimethamine_ic50.csv')
            else:
                # self.ic50_data = self.make_datapath_absolute(ic50_data)
                self.ic50_data = dir_manager.make_datapath_absolute(ic50_data)
            
            # load the data
            self.drugless_rates = self.load_fitness(self.drugless_data)
            self.ic50 = self.load_fitness(self.ic50_data)
            
            self.max_replication_rate = max(self.drugless_rates)
            self.n_genotype = self.drugless_rates.shape[0]
        
        # load fitness landscape from excel file
        elif fitness_data == 'manual':
            self.landscape_path = landscape_path
            self.landscape_data = self.load_fitness(self.landscape_path)
            
            self.max_replication_rate = max(self.landscape_data)
            self.n_genotype = self.landscape_data.shape[0]
            
        # Initial number of cells (default = 10,000 at 0000)
        if init_counts is None:
            self.init_counts = np.zeros(self.n_genotype)
            self.init_counts[0] = 10**4
        else:
            self.init_counts = init_counts
            
        self.n_allele = int(np.log2(self.n_genotype))
        
        if self.constant_pop:
            self.init_counts = self.init_counts*self.max_cells/sum(self.init_counts)
            self.init_counts = np.floor(self.init_counts)
            self.carrying_cap = False
        
        # Dose parameters
        self.curve_type = curve_type # linear, constant, heaviside, pharm, pulsed
        
        # Pharmacological paramters
        if k_abs < k_elim:
            raise Exception('Inappropriate pharmacokinetic values: k_abs < k_elim.')            
        
        self.k_elim = k_elim
        self.k_abs = k_abs 
        self.pad_right = pad_right
        self.max_dose = max_dose
        
        if slope is None:
            self.slope = self.max_dose/self.n_timestep # Ramped parameter (what time step to reach maximum dose, determines slope)
        else:
            self.slope = slope
        
        self.dose_schedule= dose_schedule
        self.prob_drop = prob_drop # probability of dropping a dose
        self.h_step = h_step # when to turn on heaviside function
        self.min_dose = min_dose 
        self.duty_cycle = duty_cycle
        
        # Generate drug dosage curves if one is not specified
        if drug_curve is None:
            self.drug_curve,u = self.gen_curves()
        else:
            self.drug_curve = drug_curve
        
        # self.static_landscape = static_landscape
        self.static_topology = static_topology
        self.static_topo_dose = static_topo_dose
            
        # Visualization parameters
        self.plot = plot # boolean
        self.plot_entropy = plot_entropy
        self.drug_log_scale = drug_log_scale # plot drugs on log scale
        self.counts_log_scale = counts_log_scale # plot counts on log scale
        self.fig_title = fig_title
        self.counts_log_scale = counts_log_scale
        if x_lim is None:
            self.x_lim = n_timestep
        else:
           self.x_lim = x_lim 
        self.y_lim = y_lim
        self.entropy_lim = entropy_lim
###############################################################################       
        
    # Load data
    # also use to load ic50 and drugless growth rate (anything from a csv)
    def load_fitness(self,data_path):
        fitness = pd.read_csv(data_path)
        cols = list(fitness.columns)
        fit_array = np.array(cols)
        fit_array = fit_array.astype('float')
        return fit_array
        
###############################################################################
    # ABM helper methods
    
    # converts decimals to binary
    def int_to_binary(self,num):
        pad = int(math.log(self.n_genotype,2))
        return bin(num)[2:].zfill(pad)
    
    # computes hamming distance between two genotypes
    def hammingDistance(self,s1,s2):
        assert len(s1) == len(s2)
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    
    # converts an integer to a genotype and padding to the left by 0s
    def convertIntToGenotype(self,anInt,pad):
    	offset = 2**pad
    	return [int(x) for x in bin(offset+anInt)[3:]]
    
    def random_mutations(self,N):
        trans_mat = np.zeros([N,N])
        for mm in range(N):
            for nn in range(N):
                trans_mat[mm, nn] = self.hammingDistance( self.int_to_binary(mm) , self.int_to_binary(nn))
        trans_mat[trans_mat>1] = 0
        trans_mat = trans_mat/trans_mat.sum(axis=1)
        return trans_mat
    
    # check if the most fit mutant is the most prevalent
    def check_stop_cond(self,counts,mm):
        final_landscape = self.gen_fit_land(self.max_dose)
        fittest_genotype = final_landscape.argmax()
        
        most_frequent_genotype = counts.argmax()
        stop_cond = False
        
        if fittest_genotype == most_frequent_genotype:
            stop_cond = True
        
        if mm >= self.n_timestep:
            raise Warning('Stop condition not reached. Increase n_timestep or adjust model parameters.')
            stop_cond = True
            
        return stop_cond

##############################################################################
    # core evolutionary model
    
    def abm(self,mm,n_genotype,P,counts):
        
        conc = self.drug_curve[mm]
            
        # __gen_fl_for_abm automatically considers carrying capacity, but
        # it does not consider timestep scale
        fit_land = self.__gen_fl_for_abm(conc, counts)
        
        fit_land = fit_land*self.timestep_scale
        death_rate = self.death_rate*self.timestep_scale
        mut_rate = self.mut_rate*self.timestep_scale
            
        if self.debug and np.mod(mm,10) == 0:
            print(str(mm))
            print(str(counts))
            print(str(fit_land))
            
        # if mm > 0 and self.bottleneck == True and np.mod(mm-self.duty_cycle*self.dose_schedule,self.dose_schedule) == 0:
        #     counts = np.round(0.1*counts)
         
        counts_t = counts

        # Kill cells
        
        counts_t = counts_t - np.random.poisson(counts*death_rate)
        
        # make sure there aren't negative numbers
        
        neg_indx = counts_t < 0
        counts_t[neg_indx] = 0
        
        # Divide cells

        daughter_counts = np.random.poisson(counts_t*fit_land)
        
        for genotype in np.arange(n_genotype):

            n_mut = np.random.poisson(daughter_counts[genotype]*mut_rate)
            # if n_mut > 0 and np.mod(mm,200)==0:
            #     print('timestep: ' + str(mm) + 
            #           '\ngenotype: ' + str(genotype) + 
            #           '\ndaughters: ' + str(daughter_counts[genotype]) +
            #           '\nn_mut: ' + str(n_mut))
            # Substract mutating cells from that allele
            daughter_counts[genotype] -= n_mut
                        
            mutations = np.random.choice(n_genotype, size=n_mut, p=P[:,genotype]).astype(np.uint8)

            # Add mutating cell to their final types
            counts_t +=np.bincount( mutations , minlength=n_genotype)
            
            # Substract mutating cells from that allele
            daughter_counts[genotype] -=n_mut

        counts_t += daughter_counts
        
        # Normalize to constant population            
        if self.constant_pop:
            # c = counts_t.astype('float') # prevent overflow or underflow errors
            # cur_size = np.sum(c)
            # c = c/cur_size
            # c = c*self.max_cells
            # counts_t = c.astype('int')
            scale = self.max_cells/np.sum(counts_t)
            counts_t = counts_t*scale
            counts_t = np.ceil(counts_t).astype('int')
        
        return counts_t
    
    def run_abm(self):
        
        n_genotype = self.n_genotype
        
        # Get transition matrix
        P = self.random_mutations( n_genotype )
        
        mm = 0
        
        # Two main modes:
        # Stop condition: run until the population reaches fixation
        # Default: run for n_timestep
        
        if self.stop_condition:
            counts = np.zeros( [1,n_genotype] , dtype=int)
            counts[0,:] = self.init_counts
            stop_condition = False
            
            while not stop_condition:
                counts_t = self.abm(mm,n_genotype,P,counts[mm])
                if len(counts.shape) == 1:
                    counts = np.append([counts],[counts_t],axis=0)
                else:
                    counts = np.append(counts,[counts_t],axis=0)
                    mm+=1
                stop_condition = self.check_stop_cond(counts_t,mm)    
            
        else:
            counts = np.zeros( [self.n_timestep, n_genotype] , dtype=int)
            counts[0,:] = self.init_counts
            
            while mm < self.n_timestep - 1:
                counts[mm+1] = self.abm(mm,n_genotype,P,counts[mm])
                mm+=1
                
        return counts, mm
    
    def simulate(self):
    
        counts = np.zeros([self.n_timestep,self.n_genotype])
        avg_counts = np.zeros([self.n_timestep,self.n_genotype])
        fixation_time = []
        
        # n_survive = 0
        for i in range(self.n_sims):
            
            if self.prob_drop > 0:
                self.drug_curve,u = self.gen_curves()
            
            counts, mm = self.run_abm()
            avg_counts += counts
            fixation_time.append(mm)

            if self.plot is True:
                self.plot_timecourse(counts_t = counts)
                
        self.counts = counts
        avg_counts = avg_counts/self.n_sims
        return avg_counts, fixation_time

##############################################################################
# Wrapper methods for fitness

    def gen_fitness(self,allele,conc,drugless_rate,ic50):
        fit = fitness.gen_fitness(self,allele,conc,drugless_rate,ic50)
        return fit
    
    def gen_fit_land(self,conc):
        fit_land = fitness.gen_fit_land(self,conc)
        return fit_land
    
    # Private to avoid confusion with gen_fit_land
    def __gen_fl_for_abm(self,conc,counts):
        fit_land = fitness.gen_fl_for_abm(self,conc,counts)
        return fit_land
    
###############################################################################
# Wrapper methods for generating drug concentration curves

    def pharm_eqn(self,t,k_elim=None,k_abs=None,max_dose=None):
        conc = pharm.pharm_eqn(self,t,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
        return conc
    
    def convolve_pharm(self,u):
        conv = pharm.convolve_pharm(self,u)
        return conv
    
    def gen_impulses(self):
        u = pharm.gen_impulses(self)
        return u
    
    def gen_on_off(self,duty_cycle=None):
        u = pharm.gen_on_off(self,duty_cycle=duty_cycle)
        return u
    
    def gen_curves(self):
        curve, u = pharm.gen_curves(self)
        return curve, u
##############################################################################
# Wrapper methods for plotting
  
    def plot_timecourse(self,counts_t=None,title_t=None):
        fig = plotter.plot_timecourse(self,counts_t=counts_t,title_t=title_t)
        return fig
    
    def plot_fitness_curves(self,fig_title='',plot_r0 = False,save=False,savename=None):
        fig = plotter.plot_fitness_curves(self,fig_title=fig_title,plot_r0 = plot_r0,save=save,savename=savename)
        return fig
    
    
    
# p2 = Population(max_cells=10**9,mut_rate=10**-3)
# np.random.seed(10)
# p2.simulate()
    