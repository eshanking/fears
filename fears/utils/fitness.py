import numpy as np
from scipy.optimize import curve_fit
# from fears.population import Population


def gen_fitness_curves(pop,conc=None):


    if conc is None:
        conc = np.logspace(-3,5,num=1000)
    
    n_genotype = pop.n_genotype

    fc = {}
    for g in range(n_genotype):
        f = np.zeros(len(conc))
        i = 0
        for c in conc:
            f[i] = gen_fitness(g,c,pop=pop) - pop.death_rate
            i+=1
        fc[g] = f

    return fc

# compute fitness given a drug concentration
def gen_fitness(pop,genotype,conc,
                drugless_rate=None,ic50=None):        

    if drugless_rate is None:
        drugless_rate = pop.drugless_rates
    if ic50 is None:
        ic50 = pop.ic50

    # logistic equation from Ogbunugafor 2016
    conc = conc/10**6 # concentration in uM, convert to M
    c = -.6824968 # empirical curve fit
    log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
    if conc <= 0:
        fitness = drugless_rate[genotype]
    else:
        fitness = log_eqn(drugless_rate[genotype],ic50[genotype])

    return fitness

def logistic_equation(conc,drugless_rate,ic50):
    """
    Logistic equation from ogbunugafor et al, PLOS CB, 2016

    Parameters
    ----------
    dugless_rate : float
        Drugless growth rate of genotype.
    ic50 : float
        ic50 of genotype.
    conc : float
        Drug concentration (in Molarity (M)).
    c : float, optional
        Logistic curve steepness parameter. The default is -0.6824968.

    Returns
    -------
    f : float
        Replication rate.

    """
    c=-0.6824968
    conc = conc/10**6
    f = drugless_rate/(1+np.exp((ic50-np.log10(conc))/c))
    
    return f

def gen_static_landscape(pop,conc):
    
    # get final landscape and seascape
    landscape = np.zeros(pop.n_genotype)
    for kk in range(pop.n_genotype):
        landscape[kk] = gen_fitness(pop,
                                    kk,
                                    pop.static_topo_dose,
                                    pop.drugless_rates,
                                    pop.ic50)
    
    if min(landscape) == 0:
        zero_indx_land = np.argwhere(landscape==0)
        landscape_t = np.delete(landscape,zero_indx_land)
        min_landscape = min(landscape_t)
    else:
        min_landscape = min(landscape)
        zero_indx_land = []
    
    seascape = np.zeros(pop.n_genotype)
    for gen in range(pop.n_genotype):
        seascape[gen] = gen_fitness(gen,conc,pop.drugless_rates,pop.ic50,pop=pop)
        
    if min(seascape) == 0:
        zero_indx_sea = np.argwhere(seascape==0)
        seascape_t = np.delete(seascape,zero_indx_sea)
        min_seascape = min(seascape_t)
    else:
        min_seascape = min(seascape)
        
    landscape = landscape - min_landscape
    landscape = landscape/max(landscape)
    
    rng = max(seascape) - min_seascape
    
    landscape = landscape*rng + min_seascape
    
    landscape[zero_indx_land] = 0
    return landscape

def gen_digital_seascape(pop,conc,gen,min_fitness=0):
    
    if pop.mic_estimate is not None:
        mic = est_mic(gen,Kmic=pop.mic_estimate,pop=pop)
    else:
        mic = est_mic(gen,growth_rate=pop.death_rate,pop=pop)
    
    if conc >= mic:
        fitness = min_fitness
    else:
        fitness = pop.drugless_rates[gen]
    return fitness

def gen_fit_land(pop,conc,mode=None):

    fit_land = np.zeros(pop.n_genotype)
            
    if pop.fitness_data == 'manual' or mode=='manual':
        fit_land = pop.landscape_data/pop.doubling_time

    else:
        
        if pop.landscape_type == 'static':
            fit_land = gen_static_landscape(pop,conc)
            
        if pop.landscape_type == 'digital':
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_digital_seascape(pop, conc, kk)
            
        elif pop.landscape_type == 'natural':
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_fitness(pop,
                                        kk,
                                        conc,
                                        pop.drugless_rates,
                                        pop.ic50)/pop.doubling_time
    
    return fit_land

# Generate fitness landscape for use in the abm method
def gen_fl_for_abm(pop,conc,counts):

    fit_land = gen_fit_land(pop,conc)
    
    # # takes the landscape at the max dose and scales the replication rate
    # # according to drug concentration
    # if pop.static_landscape:
    #     # max_fitness = max(fit_land)
    #     # fit_land = pop.gen_fit_land(pop.max_dose)
    #     # fit_land = fit_land*max_fitness/max(fit_land)
    #     fit_land = gen_fit_land(pop,conc)
    
    # if pop.static_topology:
    #     fit_land = gen_fit_land(pop,conc)
    
    # Scale division rates based on carrying capacity
    if pop.use_carrying_cap:
        division_scale = 1-np.sum(counts)/pop.carrying_cap
        if counts.sum()>pop.max_cells:
            division_scale = 0
    else:
        division_scale = 1
    
    fit_land = fit_land*division_scale
    
    return fit_land

def gen_random_seascape(n_allele,
                        drugless_limits=[1,1.5],
                        ic50_limits=[-6.5,-1.5]):
    
    n_genotype = 2**n_allele
    
    drugless_rates = np.random.uniform(min(drugless_limits),
                                    max(drugless_limits),
                                    n_genotype)
    
    ic50 = np.random.uniform(min(ic50_limits),
                            max(ic50_limits),
                            n_genotype)
    
    return drugless_rates,ic50

def randomize_seascape(pop,
                    drugless_limits=[1,1.5],
                    ic50_limits=[-6.5,-1.5]):

    n_genotype = pop.n_genotype
    
    pop.drugless_rates = np.random.uniform(min(drugless_limits),
                                        max(drugless_limits),
                                        n_genotype)

    pop.ic50 = np.random.uniform(min(ic50_limits),
                                max(ic50_limits),
                                n_genotype)
    
def fit_logistic_curve(xdata,ydata):
    
    popt,var = curve_fit(logistic_equation,xdata,ydata)
    
    return popt

def gen_null_seascape(pop,conc):

    landscape = gen_fit_land(conc,pop=pop)
    start_rates = gen_fit_land(10**-3,pop=pop)
    final_rates = gen_fit_land(10**5,pop=pop)
    # mid_rates = gen_fit_land(pop,10**1)
    
    start_points = scale_and_ignore_zeros(landscape,start_rates)
    end_points = scale_and_ignore_zeros(landscape,final_rates)
    # mid_points = scale_and_ignore_zeros(landscape,mid_rates)
    mid_points = landscape
    
    xdata = [10**-3,conc,10**5]
    
    ic50_new = []
    drugless_rates_new = []
    
    for genotype in range(len(landscape)):
        ydata = [start_points[genotype],
                mid_points[genotype],
                end_points[genotype]]
        params = fit_logistic_curve(xdata,ydata)
        ic50_new.append(params[1])
        drugless_rates_new.append(params[0])
    # find the null landscape drugless rates
    
    drugless_rates_new = scale_and_ignore_zeros(drugless_rates_new,
                                                pop.drugless_rates)
    
    return drugless_rates_new,ic50_new

def scale_and_ignore_zeros(data,target):
    """
    Scale data to range of target while ignoring the zero values in data and
    target.

    Parameters
    ----------
    data : numpy array
        Data to be scaled to the range of target.
    target : numpy array
        Target data range.

    Returns
    -------
    scaled_data : numpy array
        Scaled data to range of target. Zero values in data are set to zero
        in scaled_data and zero values in target are not used to calculate
        range.

    """
    # make sure inputs are numpy arrays
    
    if not isinstance(data,np.ndarray):
        data=np.array(data)
    if not isinstance(target,np.ndarray):
        target=np.array(target)
    
    if min(data) == 0:
        zero_indx_data = np.argwhere(data==0)
        data_t = np.delete(data,zero_indx_data)
        min_data = min(data_t)
    else:
        min_data = min(data)
        zero_indx_data = []
        
    if min(target) == 0:
        zero_indx_target = np.argwhere(target==0)
        target_t = np.delete(target,zero_indx_target)
        min_target = min(target_t)
    else:
        min_target = min(target)
        
    data = data - min_data
    data = data/max(data)

    rng = max(target) - min_target
    
    scaled_data = data*rng + min_target

    scaled_data[zero_indx_data] = 0
    
    return scaled_data

def est_mic(pop,gen,Kmic=None,growth_rate=None):
    """
    est_mic: estimates the mic based on a given Kmic (ratio of growth rate to 
    max growth rate at MIC) or based on a given growth rate.

    Parameters
    ----------
    pop : population class object
        
    gen : int
        Genotype under consideration.
    Kmic : float, optional
        Ratio of growth rate to max growth rate at MIC. The default is None.
    growth_rate : float, optional
        Growth rate at MIC. The default is None.

    Raises
    ------
    Exception
        Function requires Kmic OR growth_rate to calculate MIC.

    Returns
    -------
    mic : float
        MIC at a given growth rate or Kmic.

    """

    if Kmic is None:
        if growth_rate is None:
            raise Exception('Need a growth rate or Kmic threshold to estimate mic.')
        else:
            Kmic = growth_rate/pop.drugless_rates[gen]
    c=-0.6824968
    mic = 10**(pop.ic50[gen]+6 - c*np.log((1/Kmic)-1))
    return mic     

def sl_to_fitness(pop,g,conc,hc=None):
    """Seascape library to fitness value (in units per second)

    Args:
        pop (population class object): population class object
        g (int): genotype
        conc (float): drug concentration

    Returns:
        float: fitness
    """

    ic50 = pop.seascape_lib[str(g)]['ic50']
    
    g_drugless = pop.seascape_lib[str(g)]['g_drugless']

    if hc is None:
        hc = pop.seascape_lib[str(g)]['hill_coeff']
    
    # hc = -0.28

    f = logistic_pharm_curve(conc,ic50,g_drugless,hc)
    return f

def logistic_pharm_curve(x,IC50,g_drugless,hill_coeff):
    """Logistic dose-response curve. use if input is a single drug concentration

    Args:
        x (float): drug concentration scalar
        IC50 (float)): IC50
        g_drugless (float): drugless growth rate
        hill_coeff (float): Hill coefficient

    Returns:
        numpy array: array of growth rates
    """
    if x == 0:
        g = g_drugless
    else:
        g = g_drugless/(1+np.exp((IC50-np.log10(x))/hill_coeff))

    return g