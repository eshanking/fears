import numpy as np

# Methods for generating drug curves

# Equation for a simple 1 compartment pharmacokinetic model
def pharm_eqn(pop,t,k_elim=None,k_abs=None,max_dose=None):
    
    # scale constants according to timestep scale
    
    if k_elim is None:
        k_elim = pop.k_elim
    if k_abs is None:
        k_abs = pop.k_abs
    if max_dose is None:
        max_dose = pop.max_dose
    
    k_elim = k_elim*pop.timestep_scale
    k_abs = k_abs*pop.timestep_scale
    
    if k_elim == 0:
        conc = 1 - np.exp(-k_abs*t)
        conc = conc*max_dose
    else:
        conc = np.exp(-k_elim*t)-np.exp(-k_abs*t) 
        t_max = np.log(k_elim/k_abs)/(k_elim-k_abs)
        conc = conc/(np.exp(-k_elim*t_max)-np.exp(-k_abs*t_max))
        conc = conc*max_dose
    return conc

# Convolve dose regimen u with pharmacokinetic model
def convolve_pharm(pop,u):
                   # k_elim=0.01,
                   # k_abs=0.1,
                   # max_dose=1):
    k_elim = pop.k_elim
    k_abs = pop.k_abs
    max_dose = pop.max_dose
    
    pharm = np.zeros(pop.n_timestep)
    for i in range(pop.n_timestep):
        pharm[i] = pop.pharm_eqn(i,k_elim=k_elim,k_abs=k_abs,max_dose=max_dose)
    
    conv = np.convolve(u,pharm)
    conv = conv[0:pop.n_timestep]
    return conv

# Generates an impulse train to input to convolve_pharm()
def gen_impulses(pop):
    
    u = np.zeros(pop.n_timestep)
    impulse_indx = [0]
    
    i = 0
    
    # generate the drug dose regimen
    while impulse_indx[i] < pop.n_timestep*pop.timestep_scale-pop.dose_schedule:
        impulse_indx.append(pop.dose_schedule*(i+1))
        i+=1
    
    impulse_indx = np.array(impulse_indx)/pop.timestep_scale
    
    # eliminate random doses
    keep_indx = np.random.rand(len(impulse_indx)) > pop.prob_drop
    
    impulse_indx = impulse_indx[keep_indx]
    
    impulse_indx = impulse_indx.astype(int)
    
    u[impulse_indx]=1 # list of impulses at the simulated drug dosage times
    return u

# method to simulate an evolutionary on_off by pulsing drug 
# concentration
def gen_on_off_regimen(pop,duty_cycle=None):
    
    if duty_cycle is None:
        duty_cycle= pop.duty_cycle
    if duty_cycle is None:
        duty_cycle = 0.5
        
    u = np.zeros(pop.n_timestep)
    on = False
    for i in range(pop.n_timestep):
        if np.mod(i,pop.dose_schedule/pop.timestep_scale) == 0:
            on = True
            off_time = i + round((pop.dose_schedule*duty_cycle))/pop.timestep_scale
        if i == off_time:
            on = False
        if on:
            u[i] = pop.max_dose
    return u

# generates drug concentration curves
def gen_curves(pop):
    curve = np.zeros(pop.n_timestep)
    # print('hi')
    u = None
    if pop.curve_type == 'linear': # aka ramp linearly till timestep defined by steepness
        # cur_dose = 0
        for i in range(pop.n_timestep):
            # if i <= pop.steepness:
            #     slope = (pop.max_dose-10**(-3))/pop.steepness
            #     conc = slope*i+10**-3
            # else:
            #     # step = pop.steepness
            #     slope = (pop.max_dose-10**(-3))/pop.steepness
            # if cur_dose < pop.max_dose:
            conc = pop.slope*i*pop.timestep_scale
            
            if conc > pop.max_dose:
                conc=pop.max_dose
            # else:
            #     conc = pop.max_dose
            # cur_dose = conc
                # conc = slope*i+10**-3
            curve[i]=conc
            
    elif pop.curve_type == 'constant':
        curve[:] = pop.max_dose
        # print('here')

    elif pop.curve_type == 'heaviside':
        for i in range(pop.n_timestep):
            if i <= pop.h_step:
                curve[i] = pop.min_dose
            else:
                curve[i] = pop.max_dose 
    
    # Two compartment pharmacokinetic model
    elif pop.curve_type == 'pharm':
        for i in range(pop.n_timestep):
            curve[i] = pop.pharm_eqn(i)
    
    # Pulsed convolves an impulse train with the 1-compartment model
    elif pop.curve_type == 'pulsed':
        u = pop.gen_impulses()
        curve = pop.convolve_pharm(u)
        
    elif pop.curve_type == 'on_off':
        curve = pop.gen_on_off_regimen()
        
    return curve, u