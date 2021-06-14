import numpy as np

# compute fitness given a drug concentration
def gen_fitness(pop,allele,conc,drugless_rate,ic50):        
    c = -.6824968 # empirical curve fit

    # logistic equation from Ogbunugafor 2016
    conc = conc/10**6 # concentration in uM, convert to M
    
    if pop.static_topology:
        ic50_t = np.mean(ic50)
    if pop.static_landscape:
        rnge = max(drugless_rate) - min(drugless_rate)
        dr = ic50 - min(ic50)
        dr = dr/max(dr)
        dr = dr*rnge + min(drugless_rate)
        drugless_rate = dr
        ic50_t = np.mean(ic50)
    else:
        ic50_t = ic50[allele]
    
    # ic50 is already log-ed in the dataset
    log_eqn = lambda d,i: d/(1+np.exp((i-np.log10(conc))/c))
    if conc <= 0:
        fitness = drugless_rate[allele]
    else:
        fitness = log_eqn(drugless_rate[allele],ic50_t)

    return fitness

def gen_fit_land(pop,conc):
    
    fit_land = np.zeros(pop.n_genotype)
    
    if pop.fitness_data == 'generate':
        for kk in range(pop.n_genotype):
            fit_land[kk] = pop.gen_fitness(kk,conc,pop.drugless_rates,pop.ic50)/pop.doubling_time
    elif pop.fitness_data == 'manual':
        fit_land = pop.landscape_data/pop.doubling_time
    
    return fit_land

# Generate fitness landscape for use in the abm method
# Private to avoid confusion with gen_fit_land
def gen_fl_for_abm(pop,conc,counts):
    
    fit_land = pop.gen_fit_land(conc)
    
    # takes the landscape at the max dose and scales the replication rate
    # according to drug concentration
    if pop.static_landscape:
        max_fitness = max(fit_land)
        fit_land = pop.gen_fit_land(pop.max_dose)
        fit_land = fit_land*max_fitness/max(fit_land)
    
    if pop.static_topology:
        fit_land = pop.gen_fit_land(conc)
    
    # Scale division rates based on carrying capacity
    if pop.carrying_cap:
        division_scale = 1-np.sum(counts)/pop.max_cells
    else:
        division_scale = 1

    if counts.sum()>pop.max_cells:
        division_scale = 0
    
    fit_land = fit_land*division_scale
    
    return fit_land