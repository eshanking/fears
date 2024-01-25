import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# from fears.population import Population


def gen_fitness_curves(pop,conc=None):
    """Generates a dict of dose-response curves for a given population.

    Args:
        pop (population class object): Population object
        conc (array-like, optional): Array of concentrations to generate fitness for. 
        If None, generates a list of using np.logspace (-3,5,num=1000). Defaults to None.

    Returns:
        dict: dict of dose-reponse curves. 
        Each entry in the dict represents a genotype-specific dose response curve.
    """

    if conc is None:
        conc = np.logspace(-3,5,num=1000)
    
    n_genotype = pop.n_genotype

    fc = {}
    for g in range(n_genotype):
        f = np.zeros(len(conc))
        i = 0
        for c in conc:
            f[i] = gen_fitness(pop,g,c) - pop.death_rate
            i+=1
        fc[g] = f

    return fc

def pharmacodynamic_curve(c, gmax, gmin, mic, k):
    """pharmacodynamic model adapted from Foerster et al.

    Foerster, S., Unemo, M., Hathaway, L.J. et al. Time-kill curve analysis and
    pharmacodynamic modelling for in vitro evaluation of antimicrobials against Neisseria
    gonorrhoeae . BMC Microbiol 16, 216 (2016). 
    https://doi.org/10.1186/s12866-016-0838-9

    Args:
        c (float): drug concentration
        gmax (float): max growth rate
        gmin (float): min growth rate
        mic (float): estimated minimum inhibitory concentration
        k (float): hill coefficient
    """

    if type(c) == np.ndarray or type(c) == list:
        g = []
        for c_t in c:
            if c_t == 0:
                g.append(gmax)
            else:
                g.append(gmax - (((gmax-gmin)*(c_t/mic)**k)/((c_t/mic)**k-(gmin/gmax))))
        return g
    
    else:
        if c == 0:
            g = gmax
        else:
            g = gmax - (((gmax-gmin)*(c/mic)**k)/((c/mic)**k-(gmin/gmax)))
        return g

# compute fitness given a drug concentration
def gen_fitness(pop,genotype,conc,drugless_rate=None,ic50=None,hc=None,mic=None,
                death_model=None):    
    """Computes the fitness of a genotype at a given drug concentration.

    Args:
        pop (population class object): Population object containing data of interest
        genotype (int): Genotype of interest
        conc (float): drug concentration of interest
        drugless_rate (float, optional): Drugless growth rate. 
            If None, data is retrieved from population object. Defaults to None.
        ic50 (float, optional): Genotype-specific IC50. 
            If None, data is retrived from population object. Defaults to None.
        hc (float, optional): Hill coefficient.
            If None, data is retrived from population object. Defaults to None.

    Returns:
        float: Growth rate computed from pharmacodynamic equation.
    """

    # if pop.seascape_lib is not None:
    #     fitness = sl_to_fitness(pop,genotype,conc,hc=hc)
    #     # fitness = fitness*(60**2)

    if drugless_rate is None:
        drugless_rate = pop.drugless_rates
    if ic50 is None:
        ic50 = pop.ic50
    if mic is None:
        mic = pop.mic
    if death_model is None:
        death_model = pop.death_model

    if death_model == 'pharmacodynamic':

        gmax = pop.seascape_lib[str(genotype)]['gmax']
        gmin = pop.pharm_params['gmin']
        k = pop.pharm_params['k']
        mic = pop.seascape_lib[str(genotype)]['mic']

        # g = pharmacodynamic_curve(conc,gmax,gmin,mic,k)

        return pharmacodynamic_curve(conc,gmax,gmin,mic,k)
    
    else:

        # logistic equation from Ogbunugafor 2016
        # conc = conc/10**6 # concentration in uM, convert to M
        if hc is None:
            c = -.6824968 # empirical curve fit
        else:
            c = hc
        log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
        if conc <= 0:
            return drugless_rate[genotype]
        else:
            return log_eqn(drugless_rate[genotype],ic50[genotype])

    # return fitness

def logistic_equation(conc,drugless_rate,ic50,hc=-0.6824968):
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
    hc : float, optional
        Logistic curve steepness parameter. The default is -0.6824968.

    Returns
    -------
    f : float
        Replication rate.

    """
    
    # conc = conc/10**6
    f = drugless_rate/(1+np.exp((ic50-np.log10(conc))/hc))
    
    return f

def gen_static_landscape(pop,conc):
    """Generates a growth rate fitness landscape at a given drug concentration

    Args:
        pop (population class object): Population object
        conc (float): drug concentration at which to generate fitness landscape

    Returns:
        list: list of growth rates
    """
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

def gen_digital_seascape(pop,conc,gen,min_fitness=0.0):
    """Boostraps a 'digital seascape' from the continuous population seascape
       and computes the fitness at a given concentration.

    Args:
        pop (population class object): Population object
        conc (flaot): concentration at which to compute fitness
        gen (int): genotype of interest
        min_fitness (float, optional): Minimum fitness at conc higher than MIC. 
            Defaults to 0.

    Returns:
        float: growth rate (fitness)
    """
    
    if pop.mic_estimate is not None:
        mic = est_mic(gen,Kmic=pop.mic_estimate,pop=pop)
    else:
        mic = est_mic(gen,growth_rate=pop.death_rate,pop=pop)
    
    if conc >= mic:
        fitness = min_fitness
    else:
        fitness = pop.drugless_rates[gen]
    return fitness

def gen_fit_land(pop,conc,mode=None,**kwargs):
    """Generates the fitness landscape at a given drug concentration

    Args:
        pop (population class object): Population object
        conc (float): drug concentration
        counts (array-like): list of population counts

    Returns:
        list: lits of growth rates
    """

    fit_land = np.zeros(pop.n_genotype)
            
    if pop.fitness_data == 'manual' or mode=='manual':
        fit_land = pop.landscape_data/pop.growth_rate_norm

    else:
            
        if pop.digital_seascape:
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_digital_seascape(pop, conc, kk)
            
        else:
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_fitness(pop,kk,
                                           conc,**kwargs)/pop.growth_rate_norm
    
    return fit_land

# Generate fitness landscape for use in the abm method
def gen_fl_for_abm(pop,conc,counts):
    """Same as gen_fit_land but also scales by carrying capacity

    Args:
        pop (population class object): Population object
        conc (float): drug concentration
        counts (array-like): list of population counts

    Returns:
        list: lits of growth rates
    """

    fit_land = gen_fit_land(pop,conc)
    pos_indx = np.argwhere(fit_land>0)
    
    # Scale division rates based on carrying capacity
    if pop.use_carrying_cap:
        division_scale = 1-np.sum(counts)/pop.carrying_cap
        if counts.sum()>pop.carrying_cap:
            division_scale = 0
    else:
        division_scale = 1
    
    fit_land[pos_indx] = fit_land[pos_indx]*division_scale
    
    return fit_land

def gen_random_seascape(pop,
                        n_allele=None,
                        drugless_limits=None,
                        ic50_limits=None):
    """Generates random seascape data

    Args:
        pop (population class object): Population object
        n_allele (int, optional): Number of alleles. 
            If None, gets data from pop. Defaults to None.
        drugless_limits (list, optional): Range of drugless growth rates to generate random data from. 
            If None, gets data from pop. Defaults to None.
        ic50_limits (list, optional): Range of IC50 to generate random data from.
            If None, gets data from pop. Defaults to None.

    Returns:
        list: list of drugless growth rates and IC50s.
    """
    if n_allele is None:
        n_allele = pop.n_allele
    if drugless_limits is None:
        drugless_limits = pop.drugless_limits
    if ic50_limits is None:
        ic50_limits = pop.ic50_limits

    n_genotype = 2**n_allele

    drugless_rates = np.random.uniform(min(drugless_limits),
                                    max(drugless_limits),
                                    n_genotype)
    
    ic50 = np.random.uniform(min(ic50_limits),
                            max(ic50_limits),
                            n_genotype)
    
    return drugless_rates,ic50

def fit_logistic_curve(xdata,ydata):
    
    popt,var = curve_fit(logistic_equation,xdata,ydata)
    
    return popt

def gen_null_seascape(pop,conc,method='curve_fit'):
    """Generates a 'null seascape'.
       Methods:
       curve_fit: fits a new set of random curves to the range of the
       natural population seascape.
       sort: sorts IC50 values to match rank-order of pop drugless_rates

    Args:
        pop (population class object): Population object
        conc (float): drug concentration at which to retrieve null seascape rank order
        method (str, optional): Null seascape method. Defaults to 'curve_fit'.

    Returns:
        list: list of drugless growth rates and IC50s.
    """
    if method == 'curve_fit':
        if pop.fitness_data == 'estimate':
            hc = 0
            for key in pop.seascape_library.keys():
                hc += pop.seascape_library[key]['hill_coeff']
            hc = hc/(len(pop.seascape_library.keys()))

            landscape = gen_fit_land(pop,conc,hc=hc)
            start_rates = gen_fit_land(pop,10**-3,hc=hc)
            final_rates = gen_fit_land(pop,10**5,hc=hc)

        else:
            landscape = gen_fit_land(pop,conc)
            start_rates = gen_fit_land(pop,10**-3)
            final_rates = gen_fit_land(pop,10**5)
        # mid_rates = gen_fit_land(pop,10**1)
        
        # print(landscape)
        start_points = scale_and_ignore_zeros(landscape,start_rates)
        end_points = scale_and_ignore_zeros(landscape,final_rates)
        # mid_points = scale_and_ignore_zeros(landscape,mid_rates)
        mid_points = landscape
        
        xdata = [10**-3,conc,10**5]
        
        ic50_new = []
        drugless_rates_new = []

        # fig,ax = plt.subplots()
        
        for genotype in range(len(landscape)):

            ydata = [start_points[genotype],
                    mid_points[genotype],
                    end_points[genotype]]
            
            # ax.scatter(xdata,ydata,label=str(genotype))

            params = fit_logistic_curve(xdata,ydata)
            ic50_new.append(params[1])
            drugless_rates_new.append(params[0])
        # find the null landscape drugless rates
        # ax.set_xscale('log')
        # ax.legend()
        drugless_rates_new = scale_and_ignore_zeros(drugless_rates_new,
                                                    pop.drugless_rates)

        # fix the fact that genotype 3 in ogbunugafor data has zero fitness
        if hasattr(pop,'ic50_data'):
            if pop.ic50_data[-22:] == 'pyrimethamine_ic50.csv':
                ic50_new[3] = 0
                drugless_rates_new[3] = 0
    
    elif method == 'sort':
        
        landscape = gen_fit_land(pop,conc)
        
        dr =  np.array(pop.drugless_rates)
        ic50 = np.array(pop.ic50)
        
        landscape_t = landscape.argsort()
        landscape_ranks = np.empty_like(landscape_t)
        landscape_ranks[landscape_t] = np.arange(len(landscape))

        ic50_t = ic50.argsort()
        ic50_ranks = np.empty_like(ic50_t)
        ic50_ranks[ic50_t] = np.arange(len(ic50))

        dr_t = dr.argsort()
        dr_ranks = np.empty_like(dr_t)
        dr_ranks[dr_t] = np.arange(len(dr))

        ic50_new = np.zeros(len(landscape))
        drugless_rates_new = np.zeros(len(landscape))
        k = 0
        for g in landscape_ranks:
            indx = np.argwhere(ic50_ranks==g)
            indx = indx[0][0]
            ic50_new[k] = ic50[indx]

            indx = np.argwhere(dr_ranks==g)
            indx = indx[0][0]
            drugless_rates_new[k] = dr[indx]
            k+=1

        # drugless_rates_new = pop.drugless_rates

    if pop.fitness_data == 'estimate':
        i = 0
        for key in pop.seascape_lib.keys():
            pop.seascape_lib[key]['ic50'] = ic50_new[i]
            pop.seascape_lib[key]['g_drugless'] = drugless_rates_new[i]
            i+=1

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
        hc = pop.seascape_lib[str(g)]['hc']
    
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