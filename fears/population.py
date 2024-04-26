"""Simulate evolution with arbitrary fitness seascapes

Classes:
    PopParams
    Population
"""
import numpy as np
import math
from importlib_resources import files
from fears.utils import dir_manager, pharm, fitness, plotter, AutoRate
import pandas as pd

class PopParams:
    """Population parameters class

    Attributes:
        death_rate (float): death rate. Defaults to 0.1.
        mut_rate (float): mutation rate. Defaults to 10**-9.

        ic50_data_path (str): path to IC50 data. Defaults to 
            pyrimethamine_ic50.csv.
        drugless_data_path (str): path to drugless data. Defaults to 
            ogbunugafor_drugless.csv.

        plate_paths (list): list of strings to plate csv files for seascape estimation. 
            Defaults to test data.
        seascape_drug_conc (array-like): array of drug concentrations to be used to 
            compute the estimated fitness seascape
            Defaults to [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000] 
            ug/mL
        replicate_arrangement (str): arrangement of the genotype replicates on the plates. 
            Defaults to 'rows'
            'rows': each row represents a different genotype
            'columns': each column represents a different genotype
        data_cols (list): explicitely define the plate columns that include data. Expects 
            a list of lists of column names.
            Defaults to [['B','C','D','E','F'],['B','C','D','E','F'],['B','C',
                'D','E','F','G']] (for test data).
        moat (bool): if true, assumes the outer ring of wells in each plate is a moat 
            (i.e., contains no data). Defaults to True.
        hc_estimate (str): how AutoRate estimate the hill coefficient. Defaults to 
            'per_genotype'. If 'joint', estimates a single hill coefficient for the 
            whole seascape.

        drug_conc_range (array-like): defines the range of drug concentrations, mainly for 
            plotting purposes.
        
        constant_pop (bool): if true, normalizes population at each timestep to a constant 
            population size. Defaults to False.
        use_carrying_cap (bool): if true, attenuates growth rate as population approaches 
            carrying cap. Defaults to True.
        carrying_cap (int or float): Carrying capacity. Defaults to 10**10.
        init_counts (array-like): array of genotype counts to initialize simulations. 
            Defaults to 10^10 of genotype 0.
        n_allele (int): number of alleles in the model system.
        n_genotype (int): number of genotypes in the model system.

        fitness_data (str): Sets how to calculate fitness. Defaults to 'two-point'
            'two-point': program uses IC50 and drugless growth rate to parametrize 
                dose-response curve
            'estimate': estimates seascape data from plate_paths.
            'random': boostraps a random seascape
        seascape_type (str): For generating random seascapes. Defaults to 'natural'.
            'natural': no trade-off constraint. 
            'null': enforces no trade-off condition
        drugless_limits (array-like): limits of drugless growth rate for generating random 
            seascapes.
        ic50_limits (array-like): limits of ic50 for generating random seascapes.
        
        drug_unit (str): units of drug concentration for plotting purposes. Defaults to 
            '$\u03BC$M' (micro-Molar).
        fig_title (str): figure title used in plot_timecourse. Defaults to None
        plot_entropy (bool): if true, plots the population entropy over time in 
            plot_timecourse. Defaults to False.
        plot_drug_curve (bool): if true, plots the drug concentration curve over the 
            population curves in plot_timecourse. Defaults to True.
        x_lim (array-like): x limits for plot_timecourse. Defaults to none.
        y_lim (array-like): y limits for plot_timecourse. Defaults to none.
        counts_log_scale (bool): plots population counts in log scale. Defaults to False.
        drug_log_scale (bool): plots the drug concentration curve in log scale. Defaults 
            to False.

        n_timestep (int): total number of timesteps for a simulation. Defaults to 1000.
        timestep_scale (float): hours per timestep. Defaults to 1.
        passage (bool): if true, simulates passaging cells (bottleneck population every 
            passage_time).
        passage_time (int): frequency of cell passaging in units of timestep_scale.
        dilution (int): dilution factor for passaging. Defaults to 40.
        max_cells (int): if constant_pop is true, the population size is scaled to 
            max_cells every timestep.
        curve_type (str): sets the type of drug concentration curve. Defaults to 'pharm'.
            'pharm': one compartment pharmacokinetic curve
            'constant'
            'linear': linear ramp with slope set by slope attribute.
            'pulsed': simulates patient dosing
        prob_drop (float): for simulating patient nonadherence. Probability of forgetting 
            to take an individual dose. Defaults to 0.
        k_elim (float): elimination rate for 1-compartment pharmacokinetic model. Defaults 
            to 0.001.
        k_abs (float):  absorption rate for 1-compartment pharmacokinetic model. Defaults 
            to 0.01.
        pad_right (bool):
        max_dose (float): maximum drug concentration in drug concentration curve.
        dose_schedule (int): hours between doses. Defaults to 24.

        stop_condition (bool): if true, stops the simulation when the most frequent 
            genotype is also the most fit genotype.
        n_sims (int): number of times run_abm is called in simulate. Defaults to 10.
        debug (bool): if true, abm() prints some values useful for debugging.


    """

    def __init__(self,**kwargs):
        """Initializer

        Raises:
            Warning: Genotype/allele number mismatch.
        """
        self.death_rate, self.mut_rate = 0.1, 10**-9
        # self.ic50_data_path, self.drugless_data_path = 'pyrimethamine_ic50.csv','ogbunugafor_drugless.csv'

        self.seascape_path = None
        self.seascape_lib = None
  
        self.pharm_params_path = files('fears.data').joinpath('pharm_params_01172024.csv')

        p = files('fears.data').joinpath('pyrimethamine_ic50.csv')
        self.ic50_data_path = str(p)

        p = files('fears.data').joinpath('ogbunugafor_drugless.csv')
        self.drugless_data_path = str(p)

        plate_paths = ['20210929_plate1.csv','20210929_plate2.csv','20210929_plate3.csv']
        plate_paths = [files('fears.data').joinpath(p) for p in plate_paths]
        self.plate_paths = [str(p) for p in plate_paths]
        self.seascape_drug_conc = \
            [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000] #ug/mL
        self.replicate_arrangement = 'rows'
        self.data_cols = [['B','C','D','E','F'],\
            ['B','C','D','E','F'],['B','C','D','E','F','G']]

        min_dc = np.log10(self.seascape_drug_conc[1])
        max_dc = np.log10(max(self.seascape_drug_conc))
        self.drug_conc_range = [np.round(min_dc),np.round(max_dc)]

        self.death_model = 'pharmacodynamic'
        self.death_model_k = 0.644 # empirical
        self.gmin = -10**8/36.34 # cells/hr
        self.mic = None

        self.constant_pop, self.use_carrying_cap = False, True
        self.carrying_cap = 10**10
        self.growth_rate_norm = 1
        self.init_counts = None
        self.n_allele, self.n_genotype = None, None
        self.fitness_data = 'from_file' 
        self.moat = True
        self.hc_estimate = 'per_genotype'
        self.seascape_type = 'natural'
        self.drug_units = '$\mathrm{\mu}$g/mL'
        self.fig_title = None
        self.plot_entropy = False
        self.plot_drug_curve = True
        self.x_lim = None
        self.y_lim = None
        self.counts_log_scale = None
        self.drug_log_scale = False
        self.plot_pop_size = False
        
        self.n_timestep = 1000
        self.timestep_scale = 1
        self.passage = False
        self.passage_time = 24
        self.dilution = 40
        self.max_cells = 10**9

        self.curve_type = 'pharm'
        self.prob_drop = 0
        self.k_elim = 0.001
        self.k_abs = 0.01
        self.pad_right = True
        self.max_dose = 10
        self.dose_schedule = 24
        self.dwell = False
        self.dwell_time = 48
        self.regimen_length = None

        self.stop_condition = None
        self.plot = True
        self.n_sims = 10
        self.debug = False

        self.drugless_limits=[1,1.5]
        self.ic50_limits=[-3,3]

        # self.landscape_type = 'natural'
        self.digital_seascape = False

        self.null_seascape = False
        self.null_seascape_dose= 0
        self.null_seascape_method = 'curve_fit'

        self.ic50 = None
        self.drugless_rates = None


        for paramkey in self.__dict__.keys():
            for optkey in kwargs.keys():
                if paramkey == optkey:
                    td = {paramkey:kwargs.get(paramkey)}
                    self.__dict__.update(td)


class Population(PopParams):
    """Population class for simulating evolution.

    Population class inherits almost all atttributes from PopParams.

        Attributes:
        death_rate (float): death rate. Defaults to 0.1.
        mut_rate (float): mutation rate. Defaults to 10**-9.

        ic50_data_path (str): path to IC50 data. Defaults to 
            pyrimethamine_ic50.csv.
        drugless_data_path (str): path to drugless data. Defaults to 
            ogbunugafor_drugless.csv.

        plate_paths (list): list of strings to plate csv files for seascape estimation. 
            Defaults to test data.
        seascape_drug_conc (array-like): array of drug concentrations to be used to 
            compute the estimated fitness seascape
            Defaults to [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000] 
            ug/mL
        replicate_arrangement (str): arrangement of the genotype replicates on the plates. 
            Defaults to 'rows'
            'rows': each row represents a different genotype
            'columns': each column represents a different genotype
        data_cols (list): explicitely define the plate columns that include data. Expects 
            a list of lists of column names.
            Defaults to [['B','C','D','E','F'],['B','C','D','E','F'],['B','C',
                'D','E','F','G']] (for test data).
        moat (bool): if true, assumes the outer ring of wells in each plate is a moat 
            (i.e., contains no data). Defaults to True.

        drug_conc_range (array-like): defines the range of drug concentrations, mainly for 
            plotting purposes.
        
        constant_pop (bool): if true, normalizes population at each timestep to a constant 
            population size. Defaults to False.
        use_carrying_cap (bool): if true, attenuates growth rate as population approaches 
            carrying cap. Defaults to True.
        carrying_cap (int or float): Carrying capacity. Defaults to 10**10.
        init_counts (array-like): array of genotype counts to initialize simulations. 
            Defaults to 10^10 of genotype 0.
        n_allele (int): number of alleles in the model system.
        n_genotype (int): number of genotypes in the model system.

        fitness_data (str): Sets how to calculate fitness. Defaults to 'two-point'
            'two-point': program uses IC50 and drugless growth rate to parametrize 
                dose-response curve
            'estimate': estimates seascape data from plate_paths.
            'random': boostraps a random seascape
        seascape_type (str): For generating random seascapes. Defaults to 'natural'.
            'natural': no trade-off constraint. 
            'null': enforces no trade-off condition
        drugless_limits (array-like): limits of drugless growth rate for generating random 
            seascapes.
        ic50_limits (array-like): limits of ic50 for generating random seascapes.
        
        drug_unit (str): units of drug concentration for plotting purposes. Defaults to 
            '$\u03BC$M' (micro-Molar).
        fig_title (str): figure title used in plot_timecourse. Defaults to None
        plot_entropy (bool): if true, plots the population entropy over time in 
            plot_timecourse. Defaults to False.
        plot_drug_curve (bool): if true, plots the drug concentration curve over the 
            population curves in plot_timecourse. Defaults to True.
        x_lim (array-like): x limits for plot_timecourse. Defaults to none.
        y_lim (array-like): y limits for plot_timecourse. Defaults to none.
        counts_log_scale (bool): plots population counts in log scale. Defaults to False.
        drug_log_scale (bool): plots the drug concentration curve in log scale. Defaults 
            to False.

        n_timestep (int): total number of timesteps for a simulation. Defaults to 1000.
        timestep_scale (float): hours per timestep. Defaults to 1.
        passage (bool): if true, simulates passaging cells (bottleneck population every 
            passage_time).
        passage_time (int): frequency of cell passaging in units of timestep_scale.
        dilution (int): dilution factor for passaging. Defaults to 40.
        max_cells (int): if constant_pop is true, the population size is scaled to 
            max_cells every timestep.
        curve_type (str): sets the type of drug concentration curve. Defaults to 'pharm'.
            'pharm': one compartment pharmacokinetic curve
            'constant'
            'linear': linear ramp with slope set by slope attribute.
            'pulsed': simulates patient dosing
        prob_drop (float): for simulating patient nonadherence. Probability of forgetting 
            to take an individual dose. Defaults to 0.
        k_elim (float): elimination rate for 1-compartment pharmacokinetic model. Defaults 
            to 0.001.
        k_abs (float):  absorption rate for 1-compartment pharmacokinetic model. Defaults 
            to 0.01.
        pad_right (bool):
        max_dose (float): maximum drug concentration in drug concentration curve.
        dose_schedule (int): hours between doses. Defaults to 24.

        stop_condition (bool): if true, stops the simulation when the most frequent 
            genotype is also the most fit genotype.
        n_sims (int): number of times run_abm is called in simulate. Defaults to 10.
        debug (bool): if true, abm() prints some values useful for debugging.
    """

    def __init__(self,**kwargs):
        super().__init__(**kwargs)

        # initialize drug curve
        self.drug_curve = None
        self.impulses = None
        self.initialize_drug_curve()

        # initialize fitness
        self.initialize_fitness()

        # load pharmacodynamic parameters
        # self.pharm_params = pd.read_csv(self.pharm_params_path).to_dict()
        df = pd.read_csv(self.pharm_params_path)
        self.pharm_params = {}
        for key in df.keys():
            self.pharm_params[key] = df[key][0]

        if self.n_genotype is None:
            self.n_genotype = int(len(self.ic50))
        if self.n_allele is None:
            self.n_allele = int(np.log2(self.n_genotype))
        if int(self.n_allele) != int(np.log2(self.n_genotype)):
            raise Warning('Genotype/allele number mismatch')

        # initialize counts
        self.counts = np.zeros([self.n_timestep,self.n_genotype])

        if self.init_counts is None:
            self.init_counts = np.zeros(self.n_genotype)
            self.init_counts[0] = 10**6

        # initialize constant population condition
        if self.constant_pop:
            self.init_counts = \
                self.init_counts*self.max_cells/sum(self.init_counts)
            self.init_counts = np.floor(self.init_counts)
            self.use_carrying_cap = False
        
    def initialize_fitness(self):

        if self.ic50 is None and self.drugless_rates is None:

            if self.fitness_data == 'two-point':
                self.drugless_rates = \
                    dir_manager.load_fitness(self.drugless_data_path)
                self.ic50 = dir_manager.load_fitness(self.ic50_data_path)
                self.data_source = 'default data'
            
            elif self.fitness_data == 'estimate':
                
                # self.growth_rate_library = fitness.gen_growth_rate_library()
                # self.seascape_library = fitness.gen_seascape_library()

                f = str(files('fears.data').joinpath('plates'))
                e = AutoRate.Experiment(f,drug_conc=self.seascape_drug_conc,
                                        moat=self.moat,
                                        replicate_arrangement=\
                                            self.replicate_arrangement,
                                        data_cols=self.data_cols,
                                        hc_estimate=self.hc_estimate)
                e.execute()
                self.growth_rate_lib = e.growth_rate_lib
                self.seascape_lib = e.seascape_lib

                self.n_genotype = len(self.seascape_lib.keys())

                self.ic50 = np.zeros(self.n_genotype)
                self.drugless_rates = np.zeros(self.n_genotype)

                i = 0

                for key in self.seascape_lib.keys():
                    self.ic50[i] = self.seascape_lib[key]['ic50']
                    self.drugless_rates[i] = self.seascape_lib[key]['g_drugless']
                    i+=1
                    
                self.growth_rate_lib['drug_conc'] = self.seascape_drug_conc
                self.seascape_lib['drug_conc'] = self.seascape_drug_conc
                self.autorate_exp = e
            
            elif self.fitness_data == 'random':

                if self.n_allele is None:
                    self.n_allele = 2
                self.drugless_rates,self.ic50, = \
                    fitness.gen_random_seascape(self)
                self.n_genotype = 2**self.n_allele
                self.data_source = 'random'

            elif self.fitness_data == 'from_file':
                if self.seascape_path is not None:
                    if self.seascape_path.endswith('.csv'):
                        self.seascape_lib = pd.read_csv(self.seascape_path,index_col=0)
                    else:
                        self.seascape_lib = pd.read_excel(self.seascape_path,index_col=0)
                    self.data_source = self.seascape_path
                else: # default data
                    self.data_source = 'example data'
                    self.seascape_path = \
                        files('fears.data').joinpath('seascape_lib_09252023.csv')
                    self.seascape_lib = pd.read_csv(self.seascape_path,index_col=0)
                
                self.n_genotype = len(self.seascape_lib.keys())

                self.ic50 = np.zeros(self.n_genotype)
                self.drugless_rates = np.zeros(self.n_genotype)

                i = 0

                for key in self.seascape_lib.keys():
                    self.ic50[i] = self.seascape_lib[key]['ic50']
                    self.drugless_rates[i] = self.seascape_lib[key]['g_drugless']
                    i+=1 

        else:
            self.data_source = 'user-defined'

        if self.null_seascape:
            self.set_null_seascape(self.null_seascape_dose,self.null_seascape_method)
            
    def initialize_drug_curve(self):
        curve,u = pharm.gen_curves(self)
        self.drug_curve = curve
        self.impulses = u

    ###########################################################################
    # ABM helper methods
    def gen_neighbors(self,genotype):
        mut = range(self.n_allele)
        neighbors = [genotype ^ (1 << m) for m in mut]

        return neighbors
    
    # converts decimals to binary
    def int_to_binary(self,num):
        """
        Converts an integer to binary representation with the number of 
        digits equal to the number of alleles in the model.

        Parameters
        ----------
        num : int
            Number to be converted.

        Returns
        -------
        str
            Binary representation.

        """
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
                trans_mat[mm, nn] = self.hammingDistance( 
                    self.int_to_binary(mm) , self.int_to_binary(nn))

        trans_mat[trans_mat>1] = 0
        trans_mat = trans_mat/trans_mat.sum(axis=1)
        return trans_mat

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

    def passage_cells(self,mm,counts):
        """
        If self.passage is true, dilute cells according to self.dilution when
        the timestep is a multiple of self.passage_time.

        Parameters
        ----------
        mm : int
            Timestep.
        counts : numpy array
            Matrix of simulated cell counts.

        Returns
        -------
        counts : numpy array
            Matrix of simulated cell counts; diluted if the timestep is
            appropriate.

        """
        
        if (np.mod(mm*self.timestep_scale,self.passage_time) == 0 
            and not mm == 0 and self.passage):

            counts = np.divide(counts,self.dilution)
            counts[counts<1] == 0
            counts = [int(c) for c in counts]
            counts = np.array(counts)

        return counts

    ###########################################################################
    # core evolutionary model
    
    def abm(self,mm,n_genotype,P,counts):
        
        conc = self.drug_curve[mm]
            
        # gen_fl_for_abm automatically considers carrying capacity, but
        # it does not consider timestep scale

        fit_land = fitness.gen_fl_for_abm(self,conc,counts)
        
        fit_land = fit_land*self.timestep_scale
        death_rate = self.death_rate*self.timestep_scale
        mut_rate = self.mut_rate*self.timestep_scale
            
        if self.debug and np.mod(mm,10) == 0:
            print(str(mm))
            print(str(counts))
            print(str(fit_land))
         
        # Passage cells
        
        counts = self.passage_cells(mm, counts)
        
        counts_t = counts

        if self.death_model == 'pharmacodynamic':
            negative_fitness = fit_land < 0
            fit_land = np.abs(fit_land)
            delta_cells = np.random.poisson(counts_t*fit_land)
            delta_cells[negative_fitness] = -1*delta_cells[negative_fitness]
            counts_t = counts_t + delta_cells

        # else:
        counts_t = counts_t - np.random.poisson(counts*death_rate) # background turnover
        
        # Make sure there aren't negative numbers
        
        neg_indx = counts_t < 0
        counts_t[neg_indx] = 0
        
        # Divide cells

        daughter_counts = np.random.poisson(counts_t*fit_land)
        
        for genotype in np.arange(n_genotype):
            
            n_mut = np.random.poisson(
                daughter_counts[genotype]*mut_rate*self.n_allele)

            # Substract mutating cells from that allele
            daughter_counts[genotype] -= n_mut
            
            # Mutate cells
            mutations = np.random.choice(n_genotype, 
                                         size=n_mut, 
                                         p=P[:,genotype]).astype(np.uint8)

            # Add mutating cell to their final types
            counts_t += np.bincount( mutations , minlength=n_genotype )

        counts_t += daughter_counts

        # Normalize to constant population            
        if self.constant_pop:
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
            
            counts, mm = self.run_abm()
            avg_counts += counts
            fixation_time.append(mm)

            if self.plot is True:
                # print(type(counts))
                self.plot_timecourse(counts_t=counts)
        
        avg_counts = avg_counts/self.n_sims
        self.counts = avg_counts
        return avg_counts, fixation_time

    ###########################################################################
    # wrapper methods for plotting
    def plot_timecourse(self,**kwargs):
        fig = plotter.plot_timecourse(self,**kwargs)
        return fig

    def plot_fitness_curves(self,**kwargs):
        fig,ax = plotter.plot_fitness_curves(self,**kwargs)
        return fig,ax
    
    def plot_landscape(self,**kwargs):
        fig,ax = plotter.plot_landscape(self,**kwargs)
        return fig,ax

    ###########################################################################
    # wrapper methods for fitness

    # def __gen_fl_for_abm(self,conc,counts):
    #     """
    #     Return the fitness landscape apropriately scaled according to the 
    #     population size and carrying capacity

    #     Parameters
    #     ----------
    #     conc (float) : drug concentration
    #     counts (list) : vector of genotype population counts

    #     Returns
    #     ----------
    #     fit_land (list) : scaled fitness landscape
    #     """
    #     fit_land = fitness.gen_fl_for_abm(self,conc,counts)
    #     return fit_land

    def gen_fit_land(self,conc,**kwargs):

        fit_land = fitness.gen_fit_land(self,conc,**kwargs)
        
        return fit_land
    
    ###########################################################################
    # Wrapper methods for generating drug concentration curves

    def pharm_eqn(self,t,k_elim=None,k_abs=None,max_dose=None):
        conc = pharm.pharm_eqn(self,t,k_elim=k_elim,k_abs=k_abs,
                               max_dose=max_dose)
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
    
    def gen_passage_drug_protocol(self):
        drug_curve = pharm.gen_passage_drug_protocol(self)
        return drug_curve

    def set_drug_curve(self):
        """Sets the drug concentration curve for a given population
        """
        dc = self.gen_curves()
        self.drug_curve = dc[0]
        self.impulses = dc[1]

    ###########################################################################
    # Misc helper methods

    def reset_drug_conc_curve(self,**kwargs):
        """Resets the drug concentration curve. Also updates any paramters passed into kwargs.

           Useful when performing experiments with a large number of population objects. 
           Eliminates the need to repeatedly estimate fitness seascapes.
        """
        for paramkey in self.__dict__.keys():
            for optkey in kwargs.keys():
                if paramkey == optkey:
                    td = {paramkey:kwargs.get(paramkey)}
                    self.__dict__.update(td)
        
        self.set_drug_curve()
    
    def set_null_seascape(self,conc,method='curve_fit'):

        dr,ic50 = fitness.gen_null_seascape(self,conc,method=method)
        # print(ic50)
        self.drugless_rates,self.ic50 = dr,ic50

    def print_params(self):

        print('Biological parameters:',end='\n')
        print(' * Mutation rate: ',self.mut_rate,end='\n')
        print(' * Death rate: ',self.death_rate,end='\n')

        ic50 = [round(i,3) for i in self.ic50]

        print(' * IC50: ',end='\n')
        for g in range(len(ic50)):
            print('    ',g,': ',ic50[g],end='\n')

        print(' * Drugless growth rates: ',end='\n')
        dr = self.drugless_rates
        for g in range(len(dr)):
            print('    ',g,': ',dr[g],end='\n')

        print('Pharmacoligical parameters:',end='\n')
        print(' * Curve type: ',self.curve_type,end='\n')
        print(' * Max concentration: ',self.max_dose,end='\n')
        
        if self.curve_type == 'pharm' or self.curve_type == 'pulsed':
            print(' * k_elim: ',self.k_elim,end='\n')
            print(' * k_abs: ',self.k_abs,end='\n')

        print('Experimental parameters:',end='\n')
        print(' * N simulations: ',self.n_sims,end='\n')
        print(' * Use carrying capacity? ',self.use_carrying_cap,end='\n')
        if self.use_carrying_cap:
            print(' * Carrying capacity: ',self.carrying_cap,end='\n')

        print('Data information:',end='\n')
        print(' * Fitness data: ',self.fitness_data,end='\n')
        print(' * Data source: ',self.data_source,end='\n')
    
    ###########################################################################
    # Set wrapper method docs

    gen_fit_land.__doc__ = fitness.gen_fit_land.__doc__