import numpy as np

# compute fitness given a drug concentration
def gen_fitness(pop,allele,conc,drugless_rate,ic50):        


    # logistic equation from Ogbunugafor 2016
    conc = conc/10**6 # concentration in uM, convert to M
    c = -.6824968 # empirical curve fit
    log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
    if conc <= 0:
        fitness = drugless_rate[allele]
    else:
        fitness = log_eqn(drugless_rate[allele],ic50[allele])

    return fitness

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
    
    seascape = np.zeros(pop.n_genotype)
    for gen in range(pop.n_genotype):
        seascape[gen] = gen_fitness(pop,gen,conc,pop.drugless_rates,pop.ic50)
        
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

def gen_fit_land(pop,conc):
    
    fit_land = np.zeros(pop.n_genotype)
    
    if pop.fitness_data == 'generate':
        
        if pop.static_topology:
            fit_land = gen_static_landscape(pop,conc)
            
        else:
            for kk in range(pop.n_genotype):
                fit_land[kk] = gen_fitness(pop,kk,conc,pop.drugless_rates,pop.ic50)/pop.doubling_time
            
    elif pop.fitness_data == 'manual':
        fit_land = pop.landscape_data/pop.doubling_time
    
    return fit_land

# Generate fitness landscape for use in the abm method
# Private to avoid confusion with gen_fit_land
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
    if pop.carrying_cap:
        division_scale = 1-np.sum(counts)/pop.max_cells
        if counts.sum()>pop.max_cells:
            division_scale = 0
    else:
        division_scale = 1
    
    fit_land = fit_land*division_scale
    
    return fit_land




