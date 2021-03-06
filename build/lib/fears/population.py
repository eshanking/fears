import numpy as np
import math
import random
from importlib_resources import files
from fears.utils import dir_manager, pharm, fitness, plotter, AutoRate

class PopParams:
    """Population parameters class

    """

    def __init__(self,**kwargs):
        """Initializer
        
        Optional arguments: (all rates are in units per hour)
            death_rate (float): death rate. Defaults to 0.1.
            mut_rate (float): mutation rate. Defaults to 10**-9.

            ic50_data_path (str): path to IC50 data. Defaults to pyrimethamine_ic50.csv.
            drugless_data_path (str): path to drugless data. Defaults to ogbunugafor_drugless.csv.

            constant_pop (bool): if true, normalizes population at each timestep to a constant population size. Defaults to False.
            use_carrying_cap (bool): if true, attenuates growth rate as population approaches carrying cap. Defaults to True.
            carrying_cap (float): Carrying capacity. Defaults to 10**10.

            n_allele (int): number of alleles in the model system.
            n_genotype (int): number of genotypes in the model system.
            doubling_time (float): Average rate at which the model population divides.

            fitness_data (str): Sets how to calculate fitness. two-point: program uses IC50 and drugless growth rate to parametrize dose-response curve
            seascape_type (str): For generating random seascapes. natural: no trade-off constraint. null: enforces no trade-off condition
            
            drug_unit (str): units of drug concentration for plotting purposes.
            



        Raises:
            Warning: Genotype/allele number mismatch.
        """
        self.death_rate, self.mut_rate = 0.1, 10**-9
        # self.ic50_data_path, self.drugless_data_path = 'pyrimethamine_ic50.csv','ogbunugafor_drugless.csv'
  
        p = files('fears.data').joinpath('pyrimethamine_ic50.csv')
        self.ic50_data_path = str(p)

        p = files('fears.data').joinpath('ogbunugafor_drugless.csv')
        self.drugless_data_path = str(p)

        plate_paths = ['20210929_plate1.csv','20210929_plate2.csv','20210929_plate3.csv']
        plate_paths = [files('fears.data').joinpath(p) for p in plate_paths]
        self.plate_paths = [str(p) for p in plate_paths]
        self.seascape_drug_conc = [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000] #ug/mL

        # self.growth_rate_data = []
        # for plate_path in self.plate_paths:
        #     self.growth_rate_data.append(dir_manager.get_growth_rate_data(plate_path))

        self.constant_pop, self.use_carrying_cap = False, True
        self.carrying_cap = 10**10
        self.init_counts = None
        self.n_allele, self.n_genotype = None, None
        self.doubling_time = 1
        self.fitness_data = 'two-point' 
        self.seascape_type = 'natural'
        self.drug_units = '$\u03BC$M'
        self.fig_title = None
        self.plot_entropy = False
        self.plot_drug_curve = True
        self.x_lim = None
        self.y_lim = None
        self.counts_log_scale = None
        self.drug_log_scale = False
        
        self.n_timestep = 1000
        self.timestep_scale = 1
        self.passage = False
        self.passage_time = 24
        self.max_cells = 10**9

        self.curve_type = 'pharm'
        self.prob_drop = 0
        self.k_elim = 0.001
        self.k_abs = 0.01
        self.pad_right = True
        self.max_dose = 10
        self.dose_schedule = 24
        self.p_forget = 0

        self.stop_condition = None
        self.state = {}
        self.plot = True
        self.n_sims = 10
        self.debug = False

        self.landscape_type = 'natural'


        for paramkey in self.__dict__.keys():
            for optkey in kwargs.keys():
                if paramkey == optkey:
                    td = {paramkey:kwargs.get(paramkey)}
                    self.__dict__.update(td)
        
        # self.ic50 = dir_manager.load_fitness(self.ic50_data_path)
        # self.drugless_rates = dir_manager.load_fitness(self.drugless_data_path)
        # print(self.ic50)
    #     dr, ic50 = self.initialize_fitness()
    #     self.drugless_rates = dr
    #     self.ic50 = ic50
        
    #     if self.n_genotype is None:
    #         self.n_genotype = int(len(self.ic50))
    #     if self.n_allele is None:
    #         self.n_allele = int(np.log2(self.n_genotype))
    #     if int(self.n_allele) != int(np.log2(self.n_genotype)):
    #         raise Warning('Genotype/allele number mismatch')
        
    #     self.init_counts = np.zeros(self.n_genotype)
    #     self.init_counts[0] = 10**6
    
    # def initialize_fitness(self):
    #     if self.fitness_data == 'two-point':
    #         drugless_rates = dir_manager.load_fitness(self.drugless_data_path)
    #         ic50 = dir_manager.load_fitness(self.ic50_data_path)
    #     else:
    #         drugless_rates = None
    #         ic50 = None
    #     return drugless_rates, ic50


class Population(PopParams):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)

        # initialize constant population condition
        if self.constant_pop:
            self.init_counts = self.init_counts*self.max_cells/sum(self.init_counts)
            self.init_counts = np.floor(self.init_counts)
            self.carrying_cap = False

        # initialize drug curve
        self.drug_curve = None
        self.impulses = None
        self.initialize_drug_curve()

        # initialize fitness
        self.initialize_fitness()

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
        
    def initialize_fitness(self):
        if self.fitness_data == 'two-point':
            self.drugless_rates = dir_manager.load_fitness(self.drugless_data_path)
            self.ic50 = dir_manager.load_fitness(self.ic50_data_path)
        elif self.fitness_data == 'estimate':
            
            # self.growth_rate_library = fitness.gen_growth_rate_library()
            # self.seascape_library = fitness.gen_seascape_library()

            f = str(files('fears.data').joinpath('plates'))
            e = AutoRate.Experiment(f,drug_conc=self.seascape_drug_conc,moat=True)
            e.execute()
            self.growth_rate_lib = e.growth_rate_lib
            self.seascape_lib = e.seascape_lib

            self.ic50 = np.zeros(self.n_genotype)
            self.drugless_rates = np.zeros(self.n_genotype)

            i = 0
            print(self.seascape_lib.keys())

            for key in self.seascape_lib.keys():
                self.ic50[i] = self.seascape_lib[key]['ic50']
                self.drugless_rates[i] = self.seascape_lib[key]['g_drugless']
                i+=1
            
    def initialize_drug_curve(self):
        curve,u = pharm.gen_curves(self)
        self.drug_curve = curve
        self.impulses = u

    ###############################################################################
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
                trans_mat[mm, nn] = self.hammingDistance( self.int_to_binary(mm) , self.int_to_binary(nn))
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
        return counts

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
         
        # Passage cells
        
        counts = self.passage_cells(mm, counts)
        
        counts_t = counts

        # Kill cells
        # print(str(mm))

        counts_t = counts_t - np.random.poisson(counts*death_rate)
        
        # Make sure there aren't negative numbers
        
        neg_indx = counts_t < 0
        counts_t[neg_indx] = 0
        
        # Divide cells

        daughter_counts = np.random.poisson(counts_t*fit_land)
        
        for genotype in np.arange(n_genotype):
            
            n_mut = np.random.poisson(daughter_counts[genotype]*mut_rate*self.n_allele)

            # Substract mutating cells from that allele
            daughter_counts[genotype] -= n_mut
            
            # Mutate cells
            mutations = np.random.choice(n_genotype, size=n_mut, p=P[:,genotype]).astype(np.uint8)

            # Add mutating cell to their final types
            counts_t += np.bincount( mutations , minlength=n_genotype )

        counts_t += daughter_counts

        # Normalize to constant population            
        if self.constant_pop:
            scale = self.max_cells/np.sum(counts_t)
            counts_t = counts_t*scale
            counts_t = np.ceil(counts_t).astype('int')
        
        self.state['counts'] = counts_t
        self.state['n_mut'] = n_mut
        self.state['t'] = mm
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

    ##############################################################################
    # wrapper methods for plotting
    def plot_timecourse(self,**kwargs):
        fig,ax = plotter.plot_timecourse(self,**kwargs)
        return fig,ax

    def plot_fitness_curves(self,**kwargs):
        fig,ax = plotter.plot_fitness_curves(self,**kwargs)
        return fig,ax
    
    def plot_landscape(self,**kwargs):
        fig,ax = plotter.plot_landscape(self,**kwargs)
        return fig,ax

    ##############################################################################
    # wrapper methods for fitness

    def __gen_fl_for_abm(self,conc,counts):
        """
        Return the fitness landscape apropriately scaled according to the 
        population size and carrying capacity

        Parameters
        ----------
        conc (float) : drug concentration
        counts (list) : vector of genotype population counts

        Returns
        ----------
        fit_land (list) : scaled fitness landscape
        """
        fit_land = fitness.gen_fl_for_abm(self,conc,counts)
        return fit_land

    def gen_fit_land(self,conc,**kwargs):
        """
        Returns the fitness landscape at a given drug concentration

        Parameters
        ----------
        conc (float) : drug concentration

        Returns
        ----------
        fit_land (list) : fitness landscape

        """
        fit_land = fitness.gen_fit_land(self,conc,**kwargs)
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
    
    def gen_passage_drug_protocol(self):
        drug_curve = pharm.gen_passage_drug_protocol(self)
        return drug_curve

    def set_drug_curve(self):
        dc = self.gen_curves()
        self.drug_curve = dc[0]