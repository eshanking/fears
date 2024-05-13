import numpy as np

# Methods for generating drug curves

# Equation for a simple 1 compartment pharmacokinetic model
def pharm_eqn(pop,t,k_elim=None,k_abs=None,max_dose=None):
    """One-compartment pharmacokinetic model

    Args:
        pop (population): Population class object
        t (float or int): time
        k_elim (float, optional): Elimination rate constant. If None,
        gets data from population object. Defaults to None.
        k_abs (float, optional): Absorption rate constant. If None,
        gets data from population object. Defaults to None.
        max_dose (float, optional): Max serum drug concentration. If None,
        gets data from population object. Defaults to None.

    Returns:
        float: drug concentration according to pharmacokinetic model
    """
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
    """Models serum drug concentration of a patient undergoing a drug regimen.

    Convolves an impulse train (u) with a 1-compartment pharmacokinetic model.

    Args:
        pop (population): Population class object
        u (array-like): Impulse train of length pop.n_timestep. 0 for each time step where
        no drug is administered, 1 when drug is administered.

    Returns:
        numpy array: result of convolution
    """
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
    """Generates an impulse train of administered drug doses.

        0 for each time step whereno drug is administered, 1 when drug is administered.   
    
    Args:
        pop (population): Population class object

    Returns:
        numpy array: impulse train
    """
    u = np.zeros(pop.n_timestep)
    impulse_indx = [0]
    
    i = 0
    
    # generate the drug dose regimen

    if pop.dwell:
        dwell_time = int(pop.dwell_time/pop.timestep_scale)
    else:
        dwell_time = 0

    if pop.regimen_length is None:
        regimen_length = pop.n_timestep
    else:
        regimen_length = pop.regimen_length + pop.dwell_time

    if regimen_length > pop.n_timestep:
        regimen_length = pop.n_timestep

    while impulse_indx[i] < regimen_length*pop.timestep_scale-pop.dose_schedule:
        impulse_indx.append(pop.dose_schedule*(i+1))
        i+=1
    
    impulse_indx = np.array(impulse_indx)/pop.timestep_scale
    
    # eliminate random doses
    keep_indx = np.random.rand(len(impulse_indx)) > pop.prob_drop
    
    impulse_indx = impulse_indx[keep_indx]
    
    impulse_indx = impulse_indx.astype(int)
    
    u[impulse_indx]=1 # list of impulses at the simulated drug dosage times

    if pop.dwell:
        u[0:dwell_time] = 0

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
    """General method for generating drug concentration curves for populations

    Generates drug concentration curves based on parameters in population object

    Args:
        pop (population): Population class object

    Returns:
        list of numpy arrays: Drug concentration curve and impulse train
    """
    if pop.dwell:
        dwell_indx = int(pop.dwell_time/pop.timestep_scale)
    else:
        dwell_indx = 0
    curve = np.zeros(pop.n_timestep)
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
        curve[0:dwell_indx] = 0

    elif pop.curve_type == 'heaviside':
        for i in range(pop.n_timestep):
            if i <= pop.h_step:
                curve[i] = pop.min_dose
            else:
                curve[i] = pop.max_dose 
    
    # Two compartment pharmacokinetic model
    elif pop.curve_type == 'pharm':
        if pop.passage:
            curve = gen_passage_drug_protocol(pop)
        else:
            for i in range(pop.n_timestep):
                curve[i] = pop.pharm_eqn(i)
        # pad the curve with zeros to account for dwell time
        curve = np.pad(curve,(dwell_indx,0),'constant')
        
    
    # Pulsed convolves an impulse train with the 1-compartment model
    elif pop.curve_type == 'pulsed':
        u = pop.gen_impulses()
        curve = pop.convolve_pharm(u)
        
    elif pop.curve_type == 'on_off':
        curve = pop.gen_on_off_regimen()
        
    return curve, u

def gen_passage_drug_protocol(pop):
    """Generated drug dose over time when simulating cell passaging

    Drug concentration is constant between cell transfers (determined by pop.passage_time)

    Args:
        pop (population): Population class object

    Returns:
        numpy array: drug concentration curve
    """
    drug_curve = np.zeros(pop.n_timestep)
    
    gt = 0 # time in growth phase
    tc = 0 # time for calculating drug concentration
    
    for t in range(pop.n_timestep):
        if gt > pop.passage_time/pop.timestep_scale:
            gt = 0
            tc = t
        gt += 1
        drug_curve[t] = pharm_eqn(pop,tc)
    
    return drug_curve

def __pharm_eqn__(t,cmax,ke,ka):
    return cmax * (np.exp(-ke*t) - np.exp(-ka*t))

def est_pharm_params(thalf,tmax,t=None,ka_max=10):
    # given the half life and time of maximum concentration, estimate the parameters

    if t is None:
        t = np.linspace(0, 10, 10000)

    # estimate ke from half life
    ke = -np.log(0.5)/thalf

    # estimate ka from tmax

    ka_est = np.linspace(ke, ka_max, 100)

    est = []
    for k in ka_est:
        c = __pharm_eqn__(t, 1, ke, k)

        tmax_t = t[np.argmax(c)]

        est.append((tmax_t - tmax)**2)

    ka = ka_est[np.argmin(est)]

    return ka,ke