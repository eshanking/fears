import numpy as np
import math
import random
from fears.utils import dir_manager, pharm, fitness, plotter

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
        self.ic50_data_path, self.drugless_data_path = 'pyrimethamine_ic50.csv','ogbunugafor_drugless.csv'
        self.constant_pop, self.use_carrying_cap = False, True
        self.carrying_cap = 10**10
        self.n_allele, self.n_genotype = None, None
        self.doubling_time = 1
        self.fitness_data = 'two-point' 
        self.seascape_type = 'natural'
        self.drug_units = '$\u03BC$M'

        # load data
        self.drugless_data_path = dir_manager.make_datapath_absolute(self.drugless_data_path)
        self.ic50_data_path = dir_manager.make_datapath_absolute(self.ic50_data_path)
        
        self.init_counts = np.zeros(self.n_genotype)
        self.init_counts[0] = 10**6

        self.curve_type = 'pharm'
        self.k_elim = 0.001
        self.k_abs = 0.01
        self.pad_right = True
        self.max_dose = 10
        self.dose_schedule = 24
        self.p_forget = 0

        for paramkey in self.__dict__.keys():
            for optkey in kwargs.keys():
                if paramkey == optkey:
                    td = {paramkey:kwargs.get(paramkey)}
                    self.__dict__.update(td)
        
        
        if self.n_genotype is None:
            self.n_genotype = int(len(self.ic50))
        if self.n_allele is None:
            self.n_allele = int(np.log2(self.n_genotype))
        if int(self.n_allele) != int(np.log2(self.n_genotype)):
            raise Warning('Genotype/allele number mismatch')
        
        self.init_counts = np.zeros(self.n_genotype)
        self.init_counts[0] = 10**6


class Population(PopParams):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)

        # initialize fitness data
        self.drugless_rates = None
        self.ic50 = None
        self.initialize_fitness()

        # initialize constant population condition
        if self.constant_pop:
            self.init_counts = self.init_counts*self.max_cells/sum(self.init_counts)
            self.init_counts = np.floor(self.init_counts)
            self.carrying_cap = False

        # initialize drug curve
        self.drug_curve = None
        self.impulses = None
        self.initialize_drug_curve()
        
    def initialize_fitness(self):
        if self.fitness_data == 'two_point':
            self.drugless_rates = dir_manager.load_fitness(self.drugless_data_path)
            self.ic50 = dir_manager.load_fitness(self.ic50_data_path)
            
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
            
            # if self.prob_drop > 0:
            #     self.drug_curve,u = self.gen_curves()
            
            counts, mm = self.run_abm()
            avg_counts += counts
            fixation_time.append(mm)

            if self.plot is True:
                self.plot_timecourse(counts_t = counts)
        
        avg_counts = avg_counts/self.n_sims
        self.counts = avg_counts
        return avg_counts, fixation_time
    

# p = Population()    