import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import scipy as sp
import math
from fears.src import utils

class Population:
###############################################################################    
    # Initializer
    def __init__(self,
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
                 drug_regimen = None, # for modeling drug regimens
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
                 timestep_scale = 1,
                 x_lim = None, # plotting
                 y_lim = None # plotting
                 ):
                
        # Evolutionary parameters
        
        # Number of generations (time steps)
        self.n_timestep = n_timestep
        self.max_cells = max_cells
        
        # ABM parameters
        self.mut_rate = mut_rate
        self.death_rate = death_rate
        self.doubling_time = doubling_time
        self.timestep_scale = timestep_scale # timestep_scale = 2 -> timestep = 2 hrs, etc
        
        # Timecouse (set after running self.simulate)
        self.counts = np.zeros([self.n_timestep,16])
        self.counts_extinct = np.zeros([self.n_timestep,16])
        self.counts_survive = np.zeros([self.n_timestep,16])
                                            
        # Model parameters
        
        self.carrying_cap = True
        # self.div_scale = div_scale # Scale the division rate to simulate different organisms
        self.n_sims = n_sims # number of simulations to average together in self.simulate
        self.constant_pop = constant_pop
        # self.v2 = v2
        self.debug = debug
        
        self.fitness_data = fitness_data
        
        # Generate fitness data from IC50 and drugless growth rate data
        if fitness_data == 'generate':
            # Data paths
            if drugless_data is None:
                # self.drugless_data = "C:\\Users\\Eshan\\Documents\\python scripts\\theory division\\abm_variable_fitness\\data\\ogbunugafor_drugless.csv"
                # self.drugless_data = self.make_datapath_absolute('ogbunugafor_drugless.csv')
                self.drugless_data = utils.make_datapath_absolute('ogbunugafor_drugless.csv')
            else:
                # self.drugless_data = self.make_datapath_absolute(drugless_data)
                self.drugless_data = utils.make_datapath_absolute(drugless_data)
                
            if ic50_data is None:
                # self.ic50_data = "C:\\Users\\Eshan\\Documents\\python scripts\\theory division\\abm_variable_fitness\\data\\pyrimethamine_ic50.csv"
                # self.ic50_data = self.make_datapath_absolute('pyrimethamine_ic50.csv')
                self.ic50_data = utils.make_datapath_absolute('pyrimethamine_ic50.csv')
            else:
                # self.ic50_data = self.make_datapath_absolute(ic50_data)
                self.ic50_data = utils.make_datapath_absolute(ic50_data)
            
            # load the data
            self.drugless_rates = self.load_fitness(self.drugless_data)
            self.ic50 = self.load_fitness(self.ic50_data)
            
            self.max_replication_rate = max(self.drugless_rates)
            # determine number of alleles from data (not yet implemented)
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
            
        if drug_regimen is None:
            self.drug_regimen = u
        else:
            self.drug_regimen = drug_regimen
        
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
    
    # def make_datapath_absolute(self,filename):
    #     # takes a data file name and turns it into an absolute path
    #     p = '..' + os.sep + 'data' + os.sep + filename
    #     return p
    
    # Load data
    def load_fitness(self,data_path):
        # also use to load ic50 and drugless growth rate
        fitness = pd.read_csv(data_path)
        cols = list(fitness.columns)
        fit_array = np.array(cols)
        fit_array = fit_array.astype(np.float)
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
        # Don't include mutant no. 4
        
        trans_mat[trans_mat>1] = 0
        trans_mat = trans_mat/trans_mat.sum(axis=1)
        
    #    trans_mat[3,:] = 0
    #    trans_mat[:,3] = 0
    #    print(str(trans_mat))
        return trans_mat
    
    # compute fitness given a drug concentration
    def gen_fitness(self,allele,conc,drugless_rate,ic50):        
        c = -.6824968 # empirical curve fit

        # logistic equation from Ogbunugafor 2016
        conc = conc/10**6 # concentration in uM, convert to M
        
        # ic50 is already log-ed in the dataset
        log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
        if conc <= 0:
            fitness = drugless_rate[allele]
        else:
            fitness = log_eqn(drugless_rate[allele],ic50[allele])

        return fitness
    
    def gen_fit_land(self,conc):
        
        fit_land = np.zeros(self.n_genotype)
        
        for allele in range(self.n_genotype):
            fit_land[allele] = self.gen_fitness(allele,conc,self.drugless_rates,self.ic50)
        
        return fit_land
    
###############################################################################
    # Methods for generating drug curves
    
    # Equation for a simple 1 compartment pharmacokinetic model
    def pharm_eqn(self,t,k_elim=None,k_abs=None,max_dose=None):
        
        # scale constants according to timestep scale
        
        if k_elim is None:
            k_elim = self.k_elim
        if k_abs is None:
            k_abs = self.k_abs
        if max_dose is None:
            max_dose = self.max_dose
        
        k_elim = k_elim*self.timestep_scale
        k_abs = k_abs*self.timestep_scale
        
        conc = np.exp(-k_elim*t)-np.exp(-k_abs*t)
        t_max = np.log(k_elim/k_abs)/(k_elim-k_abs)
        conc = conc/(np.exp(-k_elim*t_max)-np.exp(-k_abs*t_max))
        conc = conc*max_dose
        return conc

    # New convolution method (much faster with numpy)
    def convolve_pharm(self,u):
                       # k_elim=0.01,
                       # k_abs=0.1,
                       # max_dose=1):
        k_elim = self.k_elim
        k_abs = self.k_abs
        max_dose = self.max_dose
        
        pharm = np.zeros(self.n_timestep)
        for i in range(self.n_timestep):
            pharm[i] = self.pharm_eqn(i,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
        
        # using FFT turns out to be much faster for convolutions!
        conv = np.convolve(u,pharm)
        conv = conv[0:self.n_timestep]
        return conv
    
    # Generates an impulse train to input to convolve_pharm()
    def gen_impulses(self):
        
        u = np.zeros(self.n_timestep)
        impulse_indx = [0]
        
        i = 0
        
        # generate the drug dose regimen
        while impulse_indx[i] < self.n_timestep*self.timestep_scale-self.dose_schedule:
            impulse_indx.append(self.dose_schedule*(i+1))
            i+=1
        
        impulse_indx = np.array(impulse_indx)/self.timestep_scale
        
        # eliminate random doses
        keep_indx = np.random.rand(len(impulse_indx)) > self.prob_drop
        
        impulse_indx = impulse_indx[keep_indx]
        
        impulse_indx = impulse_indx.astype(int)
        
        u[impulse_indx]=1 # list of impulses at the simulated drug dosage times
        return u
    
    def gen_bottleneck_regimen(self,duty_cycle=None):
        
        if duty_cycle is None:
            duty_cycle= self.duty_cycle
        if duty_cycle is None:
            duty_cycle = 0.5
            
        u = np.zeros(self.n_timestep)
        on = False
        for i in range(self.n_timestep):
            if np.mod(i,self.dose_schedule/self.timestep_scale) == 0:
                on = True
                off_time = i + round((self.dose_schedule*duty_cycle))/self.timestep_scale
            if i == off_time:
                on = False
            if on:
                u[i] = self.max_dose
        return u
    
    # generates drug concentration curves
    def gen_curves(self):
        curve = np.zeros(self.n_timestep)
        # print('hi')
        u = None
        if self.curve_type == 'linear': # aka ramp linearly till timestep defined by steepness
            # cur_dose = 0
            for i in range(self.n_timestep):
                # if i <= self.steepness:
                #     slope = (self.max_dose-10**(-3))/self.steepness
                #     conc = slope*i+10**-3
                # else:
                #     # step = self.steepness
                #     slope = (self.max_dose-10**(-3))/self.steepness
                # if cur_dose < self.max_dose:
                conc = self.slope*i*self.timestep_scale
                
                if conc > self.max_dose:
                    conc=self.max_dose
                # else:
                #     conc = self.max_dose
                # cur_dose = conc
                    # conc = slope*i+10**-3
                curve[i]=conc
                
        elif self.curve_type == 'constant':
            curve[:] = self.max_dose
            # print('here')

        elif self.curve_type == 'heaviside':
            for i in range(self.n_timestep):
                if i <= self.h_step:
                    curve[i] = self.min_dose
                else:
                    curve[i] = self.max_dose 
        
        # Two compartment pharmacokinetic model
        elif self.curve_type == 'pharm':
            for i in range(self.n_timestep):
                curve[i] = self.pharm_eqn(i)
        
        # Pulsed convolves an impulse train with the 1-compartment model
        elif self.curve_type == 'pulsed':
            u = self.gen_impulses()
            curve = self.convolve_pharm(u)
            
        elif self.curve_type == 'bottleneck':
            curve = self.gen_bottleneck_regimen()
            
        return curve, u
    
    def run_abm_v2(self):
        
        n_genotype = self.n_genotype

        # Obtain transition matrix for mutations
        P = self.random_mutations( n_genotype )

        # Keeps track of cell counts at each generation
        counts = np.zeros([self.n_timestep, n_genotype], dtype=int)

        counts[0,:] = self.init_counts
    
        for mm in range(self.n_timestep-1):
            
            if self.debug:
                if np.mod(mm,10) == 0:
                    print(str(mm))
                            
            conc = self.drug_curve[mm]
            
            fit_land = np.zeros(self.n_genotype)
            
            if self.fitness_data == 'generate':
                for kk in range(self.n_genotype):
                    fit_land[kk] = self.gen_fitness(kk,conc,self.drugless_rates,self.ic50)/self.doubling_time
                            
            elif self.fitness_data == 'manual':
                fit_land = self.landscape_data/self.doubling_time
                    
            fit_land = fit_land*self.timestep_scale
    
            # Scale division rates based on carrying capacity
            if self.carrying_cap:
                # division_scale = 1 / (1+(2*np.sum(counts[mm])/self.max_cells)**4)
                division_scale = 1-np.sum(counts[mm])/self.max_cells
            else:
                division_scale = 1
    
            if counts[mm].sum()>self.max_cells:
                division_scale = 0
            
            fit_land = fit_land*division_scale
            
            death_rate = self.death_rate*self.timestep_scale
            mut_rate = self.mut_rate*self.timestep_scale
            
            counts[mm+1] = counts[mm]
    
            # Kill cells
            
            counts[mm+1] = counts[mm+1] - np.random.poisson(counts[mm]*death_rate)
            
            # make sure there aren't negative numbers
            
            neg_indx = counts[mm+1] < 0
            counts[mm+1,neg_indx] = 0
            
            # Divide cells

            daughter_counts = np.random.poisson(counts[mm+1]*fit_land)
            
            for genotype in np.arange(n_genotype):
                # n_mut = np.random.binomial(daughter_counts[allele],mut_rate)
                n_mut = np.random.poisson(daughter_counts[genotype]*mut_rate)
                
                # Substract mutating cells from that allele
                daughter_counts[genotype] -=n_mut
                            
                mutations = np.random.choice(n_genotype, size=n_mut, p=P[:,genotype]).astype(np.uint8)
    
                # Add mutating cell to their final types
                counts[mm+1] +=np.bincount( mutations , minlength=n_genotype)
                counts[:,3] =  0
                # Substract mutating cells from that allele
                daughter_counts[genotype] -=n_mut
    
            counts[mm+1] += daughter_counts
            
            # Normalize to constant population            
            if self.constant_pop:
                cur_size = np.sum(counts[mm+1])
                counts[mm+1] = counts[mm+1]*self.init_counts[0]/cur_size
                counts[mm+1] = np.floor(counts[mm+1])

        return counts

    # Runs abm simulation n_sim times and averages results. Then sets self.counts to final result. Also quantifies survival number
    def simulate(self):
        
        counts_t = np.zeros([self.n_timestep,self.n_genotype])
        counts = np.zeros([self.n_timestep,self.n_genotype])
        counts_survive = np.zeros([self.n_timestep,self.n_genotype])
        counts_extinct = np.zeros([self.n_timestep,self.n_genotype])
        
        n_survive = 0
        for i in range(self.n_sims):
            
            if self.prob_drop > 0:
                self.drug_curve,u = self.gen_curves()
            
            counts_t = self.run_abm_v2()
                
            if any(counts_t[self.n_timestep-1,:]>0.1*self.max_cells):
                n_survive+=1
                counts_survive += counts_t
                if self.plot is True:
                    # title_t = 'Dose = ' + str(self.max_dose) + ' uM, survived'
                    self.plot_timecourse(counts_t = counts_t)
            else:
                counts_extinct += counts_t
                if self.plot is True:
                    # title_t = 'Dose = ' + str(self.max_dose) + ' uM, extinct'
                    self.plot_timecourse(counts_t = counts_t)           

            counts+=counts_t
                                                                                                      
        counts = counts/self.n_sims
        counts_survive = counts_survive/n_survive
        if  (self.n_sims - n_survive) > 0:
            counts_extinct = counts_extinct/(self.n_sims-n_survive)
        
        self.counts = counts
        self.counts_survive = counts_survive
        self.counts_extinct = counts_extinct
        
        return counts, n_survive
    
    def plot_timecourse(self,counts_t=None,title_t=None):
        
        if (self.counts == 0).all() and counts_t is None:
            print('No data to plot!')
            return
        elif counts_t is None:
            counts = self.counts
        else:
            counts = counts_t # an input other than self overrides self
        if title_t is not None:
            title = title_t
        else:
            title = self.fig_title    
            
        left = 0.1
        width = 0.8
        
        if self.plot_entropy == True:
            fig,(ax1,ax3) = plt.subplots(2,1,figsize=(6,4),sharex=True) 
            ax3.set_position([left, 0.2, width, 0.2]) # ax3 is entropy
        else:
            fig,ax1 = plt.subplots(1,1,figsize=(6,4),sharex=True)
        
        ax1.set_position([left, 0.5, width, 0.6]) # ax1 is the timecourse
                
        counts_total = np.sum(counts,axis=0)
        
        sorted_index = counts_total.argsort()
        sorted_index_big = sorted_index[-8:]
        
        colors = sns.color_palette('bright')
        colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
        
        # shuffle colors
        colors[[14,15]] = colors[[15,14]]
        
        cc = (cycler(color=colors) + 
              cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                                '--','--','--','--','--','--','--']))
        
        ax1.set_prop_cycle(cc)

        color = [0.5,0.5,0.5]
        
        if self.fitness_data == 'generate':
            ax2 = ax1.twinx() # ax2 is the drug timecourse
            ax2.set_position([left, 0.5, width, 0.6])
            ax2.set_ylabel('Drug Concentration (uM)', color=color,fontsize=20) # we already handled the x-label with ax1
            
            if self.drug_log_scale:
                if all(self.drug_curve>0):
                    drug_curve = np.log10(self.drug_curve)
                yticks = np.log10([10**-4,10**-3,10**-2,10**-1,10**0,10**1,10**2,10**3])    
                ax2.set_yticks(yticks)
                ax2.set_yticklabels(['0','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                                 '$10^1$','$10^2$','$10^3$'])
                ax2.set_ylim(-4,3)
            else:
                drug_curve = self.drug_curve
                ax2.set_ylim(0,1.1*max(drug_curve))
        
        #    ax2.plot(drug_curve, color=color, linewidth=3.0, linestyle = 'dashed')
            ax2.plot(drug_curve, color=color, linewidth=2.0)
            ax2.tick_params(axis='y', labelcolor=color)
                
            ax2.legend(['Drug Conc.'],loc=(1.25,0.93),frameon=False,fontsize=15)
            
            ax2.tick_params(labelsize=15)
        #    plt.yticks(fontsize=18)
            ax2.set_title(title,fontsize=20)
        
        # if self.normalize:
        #     counts = counts/np.max(counts)
            
        for allele in range(counts.shape[1]):
            if allele in sorted_index_big:
                ax1.plot(counts[:,allele],linewidth=3.0,label=str(self.int_to_binary(allele)))
            else:
                ax1.plot(counts[:,allele],linewidth=3.0,label=None)
                
        ax1.legend(loc=(1.25,-.12),frameon=False,fontsize=15)
    #    ax.legend(frameon=False,fontsize=15)
    #        ax.legend([str(int_to_binary(allele))])
            
        ax1.set_xlim(0,self.x_lim)
        ax1.set_facecolor(color='w')
        ax1.grid(False)
    
        # ax1.set_xlabel('Time',fontsize=20)
        ax1.set_ylabel('Cells',fontsize=20)
        ax1.tick_params(labelsize=15)
        
        if self.plot_entropy == True:
            e = self.entropy(counts)
            
            ax3.plot(e,color='black')
            ax3.set_xlabel('Time',fontsize=20)
            ax3.set_ylabel('Entropy',fontsize=20)
            if self.entropy_lim is not None:
                ax3.set_ylim(0,self.entropy_lim)
            ax3.tick_params(labelsize=15)
        
        if self.y_lim is not None:
            y_lim = self.y_lim
        else:
            y_lim = np.max(counts)
        
        if self.counts_log_scale:
            ax1.set_yscale('log')
            ax1.set_ylim(1,5*10**5)
        else:
            ax1.set_ylim(0,y_lim)
        
        xlabels = ax1.get_xticks()
        xlabels = xlabels*self.timestep_scale
        xlabels = xlabels/24
        xlabels = np.array(xlabels).astype('int')
        ax1.set_xticklabels(xlabels)
        ax1.set_xlabel('Days',fontsize=20)
            
            
        # x_labels = ax1.get_xlabels()
        # for label in xlabels:
            
        
        plt.show()
        return fig
    
    # Calculate the shannon-gibbs entropy (normalized population size)
    def entropy(self,counts=None):
        if counts is None:
            counts = self.counts

        k = np.sum(counts,1)
        entropy = np.zeros(counts.shape[0])
        counts_t = np.zeros(counts.shape)
        for i in range(counts.shape[0]):
            counts_t[i,:] = np.divide(counts[i,:],k[i])
            # counts_log = np.log(counts_t[i,:]+1)
            # entropy[i] = np.dot(counts_t[i,:],counts_log)
            entropy[i] = sp.stats.entropy(counts_t[i,:])

        return entropy
    
    def plot_fitness_curves(self,fig_title=''):
    
        drugless_rates = self.drugless_rates
        ic50 = self.ic50
        
        fig, ax = plt.subplots(figsize = (13,6))
        
        powers = np.linspace(-3,5,20)
        conc = np.power(10*np.ones(powers.shape[0]),powers)
        fit = np.zeros(conc.shape[0])
        
        colors = sns.color_palette('bright')
        colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
        colors[[14,15]] = colors[[15,14]]
        
        cc = (cycler(color=colors) + 
               cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                                '--','--','--','--','--','--','--']))
        ax.set_prop_cycle(cc) 
        
        for allele in range(16):
            if allele == 3:
                fit = np.zeros(conc.shape[0])
            if allele > 3:
                for j in range(conc.shape[0]):
                    fit[j] = self.gen_fitness(allele,conc[j],drugless_rates,ic50)
            else:
                for j in range(conc.shape[0]):
                    fit[j] = self.gen_fitness(allele,conc[j],drugless_rates,ic50)
            ax.plot(powers,fit,linewidth=3,label=str(self.int_to_binary(allele)))

        ax.legend(fontsize=15,frameon=False,loc=(1,-.10))
        ax.set_xticks([-3,-2,-1,0,1,2,3,4,5])
        ax.set_xticklabels(['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                             '$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
        
        plt.title(fig_title,fontsize=20)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        
        plt.xlabel('Drug concentration ($\mathrm{\mu}$M)',fontsize=20)
        plt.ylabel('Growth Rate',fontsize=20)
        ax.set_frame_on(False)
        
        return fig, ax
