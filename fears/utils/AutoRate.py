import pandas as pd
import os
import scipy.optimize as sciopt
import scipy
import scipy.interpolate as sciinter
import matplotlib.pyplot as plt
import numpy as np
from functools import partial
from fears.utils.plotter import gen_color_cycler
import datetime

class Experiment():
    """Experiment class for a given plate reader experiment
    """
    def __init__(self,
                folder_path,
                moat = False,
                replicate_arrangement = 'rows',
                mode = 'timeseries',
                hc_estimate = 'per_genotype',
                exp_layout_path = None,
                ref_data_path = None,
                ref_genotypes = 0,
                ref_keys = 'B2',
                t_obs=None,
                drug_conc = None,
                units = 'ug/mL',
                data_cols = None,
                debug = False):
        """Initializer
        
        Args:
            folder_path (str): path of plate reader data
            moat (bool, optional): If true, assumes the outer row of the plate is a moat. Defaults to False.
            replicate_arrangement (str, optional): Determines if replicates are arranged in rows or columns. Defaults to rows.
            drug_conc (list, optional): Drug concentrations corresponding to the drug diluation scheme. If none, defaults to 
            [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000].
            units (str, optional): Drug concentration units. Defaults to ug/mL
            data_cols (list of lists, optional): for each plate, defines which columns include genotype growth rate data
        """
        self.moat = moat
        self.folder_path = folder_path
        self.replicate_arrangement = replicate_arrangement
        self.data_cols = data_cols
        
        self.plates = []
        self.units = units
        self.debug = debug

        self.mode = mode
        self.exp_layout_path = exp_layout_path
        self.ref_data_path = ref_data_path
        self.ref_genotypes = ref_genotypes
        self.ref_keys = ref_keys
        self.t_obs = t_obs

        self.hc_estimate = hc_estimate
        
        # if mode is single_measurement, assumes each drug conc refers to a single plate
        if drug_conc is None:
            self.drug_conc = [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000]
        else:
            self.drug_conc = drug_conc

        
    def execute(self):
        """run the growth rate analyses
        """
        
        self.plate_data_paths = self.get_plate_data_paths()

        # get plate data paths
        i = 0
        for pdp in self.plate_data_paths:
            if self.data_cols is not None:
                dc = self.data_cols[i] # expects to be a list of columns
                i+=1
            else:
                dc = None
            p = Plate(pdp,self.drug_conc,
                      debug=self.debug,
                      moat=self.moat,
                      replicate_arrangement=self.replicate_arrangement,
                      data_cols=dc,
                      mode=self.mode,
                      exp_layout_path=self.exp_layout_path,
                      ref_data_path=self.ref_data_path,
                      ref_genotypes=self.ref_genotypes,
                      ref_keys=self.ref_keys,
                      t_obs=self.t_obs)
            p.execute()
            self.plates.append(p)
            
        
        # compute growth rate and seascape libraries
        self.growth_rate_lib = self.gen_growth_rate_lib()
        self.seascape_lib = self.gen_seascape_lib()

    def get_plate_data_paths(self):
        """Gets plate data paths

        Returns:
            list: list of plate data paths
        """
        plate_files = os.listdir(path=self.folder_path)

        #Need to make sure we are only attempting to load .csv or .xlsx data
        plate_files = [i for i in plate_files if ('.csv' in i) or ('.xlsx' in i)]

        plate_files.sort()

        plate_data_paths = []

        for pf in plate_files:
            if pf != '.DS_Store':
                plate_path = self.folder_path + os.sep + pf
                plate_data_paths.append(plate_path)

        return plate_data_paths
    
    def gen_growth_rate_lib(self):
        """Generates the complete growth rate library by compining the growth rate libraries of each plate

        Returns:
            dict: growth rate library
        """
        gl = {}
        if self.mode == 'timeseries':
            
            offset = 0
            for p in self.plates:

                keys = [int(k) for k in p.growth_rate_lib.keys()]
                
                for k in keys:
                    rep_num =  k + offset
                    gl[str(rep_num)] = p.growth_rate_lib[str(k)]
                
                offset += max(keys) + 1
        elif self.mode == 'single_measurement':
            # get the total number of genotypes
            gd = self.plates[0].genotype_dict
            genotypes = [g for g in gd.keys() if g.isnumeric()]
            num = 0

            for g in genotypes:
                g_dict = {}
                num = 0
                for p in self.plates:
                    key = str(self.drug_conc[num])
                    g_dict[key] = p.growth_rate_lib[g]
                    num+=1
                gl[g] = g_dict
        return gl
    
    def gen_seascape_lib(self):
        """Fits raw estimated growth rate values to a Hill dose-response curve

        Args:
            pop (population class object, optional): population class object. Defaults to None.
            debug (bool, optional): generates plots useful for debugging if true. Defaults to False.

        Raises:
            ValueError: raises error if there is no growth rate library in the experiment object

        Returns:
            dict: seascape library
        """


        if not 'growth_rate_lib' in self.__dict__:
            raise ValueError('No growth rate library in population.')
        else:
            gl = self.growth_rate_lib
        
        sl = {}
        dc = self.drug_conc

        replicates = self.growth_rate_lib.keys()

        if self.hc_estimate == 'joint':
            hc = self.estimate_hill_coeff()
        else:
            hc = None
        
        for r in replicates:

            growth_dict = gl[r]

            dc = [float(c) for c in growth_dict.keys()]
            growth_rates = [growth_dict[k] for k in growth_dict.keys()]

            popt = self.fit_hill_curve(dc,growth_rates,hc=hc)
            
            # inidices of optimized parameter vector
            ic50_indx = 0
            g_drugless_indx = 1
            hill_coeff_indx = 2

            ic50 = popt[ic50_indx]
            g_drugless = popt[g_drugless_indx]
            if hc is None:
                hill_coeff = popt[hill_coeff_indx]
            else:
                hill_coeff = hc

            d_t = {'ic50':ic50,
                'g_drugless':g_drugless,
                'hill_coeff':hill_coeff}

            sl[r] = d_t

        return sl

    def fit_hill_curve(self,xdata,ydata,hc=None):
        """Fits dose-response curve to growth rate data

        Args:
            xdata (list or numpy array): drug concentration curve from plate experiment
            ydata (list or numpy array): growth rate versus drug concetration for a given replicate

        Returns:
            list: List of optimized paramters: IC50, drugless growth rate, and Hill coefficient
        """
        # interpolate data
        xd_t = xdata
        yd_t = ydata
        f = sciinter.interp1d(xdata,ydata)

        if min(xdata) == 0:
            xmin = np.log10(xdata[1]) # if xdata starts at zero, set the new xmin to be the log of the next smallest value
        else:
            xmin = np.log10(min(xdata))
        xmax = np.log10(max(xdata))

        xdata = np.logspace(xmin,xmax)
        if not xdata[0] == 0:
            xdata = np.insert(xdata,0,0) # add zero back to xdata if removed before (because of log(0) error)

        ydata = f(xdata) # interpolate new ydata points

        if hc is None: # want to estimate hill coefficient as well
            p0 = [0,ydata[0],-0.08]

            if ydata[0] == 0:
                g_drugless_bound = [0,1]
            else:
                # want the estimated drugless growth rate to be very close to the value given in ydata
                g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

            bounds = ([-5,g_drugless_bound[0],-1],[4,g_drugless_bound[1],-0.001]) # these aren't magic numbers these are just starting parameters that happen to work

            popt, pcov = sciopt.curve_fit(self.logistic_pharm_curve_vectorized,
                                                xdata,ydata,p0=p0,bounds=bounds)
        
        else: # we already know the hill coefficient, estimate everything else
            p0 = [0,ydata[0]]

            # print(p0)

            if ydata[0] == 0:
                g_drugless_bound = [0,1]
            else:
                # want the estimated drugless growth rate to be very close to the value given in ydata
                g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

            fitfun = partial(self.logistic_pharm_curve_vectorized,hill_coeff=hc)

            bounds = ([-5,g_drugless_bound[0]],[4,g_drugless_bound[1]])
            popt, pcov = scipy.optimize.curve_fit(fitfun,
                                                xdata,ydata,p0=p0,bounds=bounds)            

        if self.debug:
            est = [self.logistic_pharm_curve(x,popt[0],popt[1],popt[2]) for x in xdata]
            fig,ax = plt.subplots()
            ax.plot(xd_t,yd_t)
            ax.plot(xdata,ydata)
            ax.plot(xdata,est)
            ax.set_xscale('log')
            ax.set_title('IC50 = ' + str(popt[0]))

        return popt
    
    def logistic_pharm_curve_vectorized(self,x,IC50,g_drugless,hill_coeff):
        """Defines the logistic dose-response curve. Use if the input is a vector of drug concentration curves

        Args:
            x (numpy array): drug concentration vector
            IC50 (float)): IC50
            g_drugless (float): drugless growth rate
            hill_coeff (float): Hill coefficient

        Returns:
            numpy array: array of growth rates
        """
        g = []

        for x_t in x:
            if x_t == 0:
                g.append(g_drugless)
            else:
                g.append(g_drugless/(1+np.exp((IC50-np.log10(x_t))/hill_coeff)))

        return g

    def save_results(self):
        seascape_df = pd.DataFrame.from_dict(self.seascape_lib)
        growth_rate_df = pd.DataFrame.from_dict(self.growth_rate_lib)

        seascape_df.to_csv(path_or_buf='sl.csv')
        growth_rate_df.to_csv(path_or_buf='gr.csv')
    
    def hill_curve_loss(self,genotype,x0,gr_lib=None,hc_est=None):
        """Calculates the loss for an individual hill curve based on the absolute value of the difference

        Args:
            genotype (int): genotype whose loss you want to calculate
            x0 (list or array): hill paramters: ic50, g_drugless, and hill coefficient (in that order)
            gr_lib (dict, optional): Growth rate library. Defaults to None.
            hc_est (float, optional): Hill coefficient estimate. Defaults to None.

        Returns:
            flaot: loss
        """
        if gr_lib is None:
            gr_lib = self.growth_rate_lib

        # drug_conc = gr_lib['drug_conc']
        drug_conc = self.drug_conc

        ic50_est = x0[0]
        g_drugless_est = x0[1]

        if hc_est is None:
            hc_est = x0[2]
        l=0
        for c in range(len(drug_conc)):

            gr_est = self.logistic_pharm_curve(float(c),ic50_est,g_drugless_est,hc_est)
            # plt.plot(gr_est)
            # gr_vect = gr_lib[str(genotype)]
            gr_vect = self.get_gr_vect_from_gr_lib(str(genotype),gr_lib=gr_lib)

            l += np.abs(gr_est - gr_vect[c])

        return l

    def get_gr_vect_from_gr_lib(self,g,gr_lib = None):

        if gr_lib is None:
            gr_lib = self.growth_rate_lib

        gr_vect_t = gr_lib[g]
        gr_vect = [gr_vect_t[key] for key in gr_vect_t.keys()]
        return gr_vect
    
    def hill_coeff_loss(self,hc):
        """Calculates the loss across the entire seascape for a given hill coefficient.

        Args:
            hc (float): Hill coefficient

        Returns:
            float: loss
        """
        # estimates the loss for different hill coefficients
        # hc = np.zeros(self.n_genotype)

        # xdata = self.growth_rate_lib['drug_conc']
        xdata = self.drug_conc
        loss = 0

        n_genotype = 0
        for key in self.growth_rate_lib.keys():
            if key.isnumeric():
                n_genotype+=1

        for g in range(n_genotype):
            # print(g)
            ydata_t = self.growth_rate_lib[str(g)]
            ydata = [ydata_t[key] for key in ydata_t.keys()]
            popt = self.fit_hill_curve(xdata,ydata,hc=hc)

            x0 = [popt[0],popt[1]]

            loss += self.hill_curve_loss(g,x0,hc_est=hc)


        # loss = self.hill_curve_loss(x0)

        return loss
    
    def estimate_hill_coeff(self):
        """Uses scipy minimize_scalar to estimate the hill coefficient for a seascape

        Returns:
            float: optimized hill coefficient
        """
        hc = scipy.optimize.minimize_scalar(self.hill_coeff_loss,
                                            bounds=[-10,-0.01],
                                            method='bounded')

        return hc.x
    
    def logistic_pharm_curve(self,x,IC50,g_drugless,hill_coeff):
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

    def plot_seascape(self):

        fig,ax = plt.subplots()

        cc = gen_color_cycler()
        ax.set_prop_cycle(cc)

        x = self.drug_conc

        for g in self.growth_rate_lib.keys():
            y = []
            yerr = []
            dose_response_curve = self.growth_rate_lib[g]
            for key in dose_response_curve.keys():
                avg = dose_response_curve[key]['avg']*3600
                err = dose_response_curve[key]['std']*3600
                y.append(avg)
                yerr.append(err)
            ax.errorbar(x,y,yerr=yerr,label=g)
        ax.set_xscale('log')
        ax.set_ylabel('Growth rate ($hr^{-1}$)')
        ax.set_xlabel('Drug concentration (ug/mL)')
        ax.legend(loc=(1.05,0),frameon=False)
        return fig


class Plate():
    """96-well plate object
    """
    def __init__(self,
                 data_path,
                 drug_conc=None,
                 mode='timeseries',
                 replicate_arrangement='rows',
                 moat=False,
                 debug=False,
                 data_cols=None,
                 ref_data_path=None,
                 exp_layout_path=None,
                 ref_genotypes='0',
                 ref_keys='B2',
                 t_obs=None,
                 tmax=None,
                 data_start = None):
        """Initializer

        Args:
            data_path (str): csv file path
            drug_conc (list of floats): drug concentration gradient
            replicates (str, optional): Determines if replicates are arranged in rows or columns. Defaults to rows.
            moat (bool, optional): If true, assumes the outer row of the plate is a moat. Defaults to False.
            debug (bool, optional): If true, plots growth curve and estimated growth curve. Defaults to False.
        """
        self.moat = moat
        self.data_path = data_path
        self.mode = mode
        self.data_cols = data_cols
        self.tmax = tmax
        self.data_start = data_start

        self.exp_layout_path = exp_layout_path
        if exp_layout_path is not None:
            self.genotype_dict = self.parse_exp_layout_file()
        else:
            self.genotype_dict = None

        if self.mode == 'timeseries':
            self.data = self.parse_data_file(data_path)

        elif self.mode == 'single_measurement':
            self.ref_genotypes = ref_genotypes
            self.ref_keys = ref_keys
            self.t_obs = t_obs
            self.data = self.parse_od_data_file(data_path)
            self.ref_data = self.parse_data_file(ref_data_path)
            self.set_background()

        if drug_conc is None:
            self.drug_conc = [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000]
        else:
            self.drug_conc = drug_conc

        self.replicate_arrangement = replicate_arrangement
        self.debug = debug
        self.n_genotype = None
        self.start_time = 0

    def execute(self):
        
        if self.mode == 'timeseries':
            if self.genotype_dict is None:
                self.background_keys = self.get_background_keys()
                self.data_keys = self.get_data_keys()
                self.growth_rate_lib = self.gen_growth_rate_lib_ts()
            else:
                self.growth_rate_lib = self.gen_growth_rate_lib_ts()
        elif self.mode == 'single_measurement':
            self.set_ref_params()
            self.set_constants()
            self.growth_rate_lib = self.gen_growth_rate_lib_sm()
            
    def set_constants(self):
        constants = self.compute_constants()
        self.OD_max = constants['OD_max']
        self.OD_0 = constants['OD_0']
        self.L = constants['L']

    def set_ref_params(self):
        self.ref_params = self.get_reference_params(self.ref_genotypes,
                                                    self.ref_keys,
                                                    self.ref_data)

    # def get_start_time(self):
         
    #     # load csv or xlsx
    #     p = self.data_path
    #     if '.csv' in p:
    #         df = pd.read_csv(p)
    #     elif '.xlsx' in p:
    #         df = pd.read_excel(p)
    #     col0 = df.keys()[0]
    #     col0 = df[col0]
    #     col0 = np.array(col0)
        
    #     if 'Time:' in col0:
            
    #         start_time_row = np.argwhere(col0=='Time:')[0][0]
    #         start_time_col = df[df.keys()[4]]
    #         start_time = start_time_col[start_time_row]

    #         if 'PM' in start_time and start_time[0:2] != '12':
    #             add_pm = 12
    #         else:
    #             add_pm = 0
            
    #         # get the semicolon
    #         semicolon_loc = start_time.index(':')
    #         start_time = 60*(int(start_time[0:semicolon_loc]) + add_pm) + \
    #             int(start_time[semicolon_loc+1:semicolon_loc+3])
    #         # print(start_time)
    #         self.start_time = start_time
        
    #     return start_time

    def get_start_time(self,col=4,df=None):

        if df is None:
            p = self.data_path
            if '.csv' in p:
                df = pd.read_csv(p)
            elif '.xlsx' in p:
                df = pd.read_excel(p)

        # first start time is shaking, so we take the second (start of scan)
        f = df[df == 'Start Time'].stack().index.tolist()[1]

        row = f[0]
        date_time = df.iloc[row,col]

        yr = int(date_time[0:4])
        mon = int(date_time[5:7])
        day = int(date_time[8:10])

        hr = int(date_time[11:13])
        min = int(date_time[14:16])
        sec = int(date_time[17:19])

        dt = datetime.datetime(yr,mon,day,hour=hr,minute=min,second=sec)


        return dt

    def parse_data_file(self,p,data_start=None):
        """Strips metadata from raw data file to obtain timeseries OD data

        Args:
            p (str): path to data file

        Returns:
            pandas dataframe: dataframe of raw data
        """
        
        # load csv or xlsx
        if '.csv' in p:
            df = pd.read_csv(p)
        elif '.xlsx' in p:
            df = pd.read_excel(p)

        # get the first column (leftmost) of the data
        # cycle nr is always in the leftmost column
        time_col = df[df.keys()[0]]
        time_array = np.array(time_col)
            
        if data_start == None:
            # raw data starts after cycle nr.
            if any(time_array == 'Cycle Nr.'):
                data_start_indx = np.argwhere(time_array == 'Cycle Nr.')
            elif any(time_array == 'Time [s]'):
                data_start_indx = np.argwhere(time_array == 'Time [s]')
            else:
                raise Exception('Unknown file format. Expected either Cycle Nr. or Time [s] as column headings.')

            #sometimes the data gets passed in very unraw
            if len(data_start_indx) == 0:
                return df
            
            data_start_indx = data_start_indx[0][0] # get scalar from numpy array
        else:

            if any(time_array == data_start):
                data_start_indx = np.argwhere(time_array == data_start)

            else:
                raise Exception('Specified data start string not found')
            
            data_start_indx = data_start_indx[0][0] + 1 # get scalar from numpy array
    

        # filter header from raw data file
        df_filt = df.loc[data_start_indx:,:]

        # change the column names
        df_filt.columns = df_filt.iloc[0] # get the top row from the dataframe and make it the column names
        df_filt = df_filt.drop(df_filt.index[0]) # remove the top row from the dataframe

        # find data end and filter empty cells
        first_col = df_filt[df_filt.keys()[0]] # get the first column of data
        x = pd.isna(first_col) # find where na occurs (empty cell)
        data_end_indx = x[x].index[0] 

        df_filt = df_filt.loc[:data_end_indx-1,:] # filter out the empty cells

        return df_filt

    def get_background_keys(self):
        """Gets the dataframe keys for the background (aka moat)

        Returns:
            list: list of background keys
        """
        # row A, row H, col 1, and col 12
        first_plate_col = 1 #first column of the plate
        last_plate_col = 12 #last column of the plate

        k = self.data.keys()

        # filter out keys that don't refer to wells
        k_filt = []
        for key in k:
            if self.check_if_key_is_well(key):
                k_filt.append(key)

        k = k_filt

        if self.moat:

            # this block of code defines the outer ring of a 96-well plate

            bg_keys = [y for y in k if int(y[1:]) == first_plate_col] # col 1
            bg_keys = bg_keys + [y for y in k if (int(y[1:]) == last_plate_col and y not in bg_keys)]
            bg_keys = bg_keys + [y for y in k if (y[0] == 'A' and y not in bg_keys)]
            bg_keys = bg_keys + [y for y in k if (y[0] == 'H' and y not in bg_keys)]
        
            if self.data_cols is not None: # add whatever is not in a data column to the background keys
                
                bgk_t = self.get_not_data_keys()
                
                bg_keys = list(set(bg_keys) | set(bgk_t)) # union avoids repeats

        elif self.data_cols is not None: # add whatever is not in a data column to the background keys
            bg_keys = self.get_not_data_keys()
        else:
            bg_keys = None

        return bg_keys
        
    def get_not_data_keys(self):
        """This function is called if data_cols is not None. Return the list of keys that don't refer to wells that have actual data.

        Returns:
            list: list of keys for wells not in data_keys
        """
        first_plate_col = 1 #first column of the plate
        last_plate_col = 12 #last column of the plate
        k = self.data.keys()

        # filter out keys that don't refer to wells
        k_filt = []
        for key in k:
            if self.check_if_key_is_well(key):
                k_filt.append(key)

        k = k_filt

        data_keys = []

        if self.replicate_arrangement == 'rows': # genotypes are defined by letters corresponding to rows
            for r in self.data_cols:
                for col in range(first_plate_col,last_plate_col):
                    key = r + str(col)
                    data_keys.append(key)
        else:
            for col in self.data_cols:
                for r in ['A','B','C','D','E','F','G','H']:
                    key = r + str(col)
                    data_keys.append(key)
        bg_keys = [key for key in k if key not in data_keys] # for every well in the data set, add it to the background key list if it is not in the data key list

        return bg_keys

    def get_data_keys(self):
        """Gets the dataframe keys for the data

        Args:
            df (pandas dataframe): datafram containing raw OD data

        Returns:
            list: list of keys
        """
        if self.background_keys is None:
            data_keys = self.data.keys()
        else:
            data_keys = [k for k in self.data.keys() if k not in self.background_keys]
        
        data_keys = [k for k in data_keys if self.check_if_key_is_well(k)]

        return data_keys

    def check_if_key_is_well(self,key):
        """Checks if key could refer to a well (i.e. key is in format 'X##' where X is a letter and # are numbers)
           For instance, returns True is key is 'A11', returns False if key is 'Time (s)' 

        Args:
            key str: key to check

        Returns:
            boolean: True if key could refer to a well, False if other.
        """
        isKey = False

        max_key_length = 3
        min_key_length = 2

        if len(key) <= max_key_length and len(key) >= min_key_length:
            if key[0].isalpha(): # is the first element a letter?
                if key[1:].isdigit(): # is the last element (or last two elements) a number?
                    isKey = True

        return isKey
    
    def est_growth_rate(self,growth_curve,t=None):
        """Estimates growth rate from OD growth curve

        Args:
            growth_curve (list or numpy array): vector of OD data
            t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

        Returns:
            float: Growth rate in units of 1/s
        """
        
        if t is None:
            t = np.arange(len(growth_curve))

        if normalize:
            norm_factor = max(growth_curve)
            growth_curve = [g/norm_factor for g in growth_curve]
        else:
            norm_factor = 1

        cc_est = np.mean(growth_curve[-2:])
        gr_est = rolling_regression(t,growth_curve)
        p0 = [10**-6,0.05,1] # starting parameters

        popt, pcov = sciopt.curve_fit(self.logistic_growth_with_lag,
                                            t,growth_curve,p0=p0,
                                            bounds=(0,1))
        
        rate_indx = 0 # index of the optimized data points referring to the growth rate
        p0_indx = 1 # index of the optimized data points referring to initial population size
        carrying_cap_indx = 2 # index of the optmized data points referring to carrying capacity

        r = popt[rate_indx]
        p0 = popt[p0_indx]
        cc = popt[carrying_cap_indx]

        min_carrying_cap = 0.1

        if r < 0: # if the growth rate is negative
            r = 0
        if cc < p0: # if the carrying capacity is less than the initial population size
            r = 0
        if cc < min_carrying_cap: # if the carrying cap is less the minimum threshold
            r = 0
        
        if self.debug:
            fig,ax = plt.subplots()

            ax.plot(t,growth_curve)

            est = self.logistic_growth_curve(t,popt[0],popt[1],popt[2])
            
            ax.plot(t,est)
            # print(popt[0])
            p0 = round(popt[1]*10**5)/10**5
            k = round(popt[2]*10**5)/10**5
            r = round(r*10**5)/10**5
            title = 'rate = ' + str(r*(60**2)) + ' cc = ' + str(k)
            ax.set_title(title)        

        return r

    def get_growth_rates_from_df(self):
        
        """Estimates the growth rates from timeseries growth data in a dataframe

        Returns:
            growth_rates: dict
                dictionary of growth rates for each experimental condition
        """

        growth_rates = {}
        df = self.data
        if self.genotype_dict is None:
            data_keys = self.get_data_keys()
            time = df['Time [s]']

            for k in data_keys:
                gr = np.array(df[k]) # growth rate time series data
                time = np.array(time) # time vector
                # fig,ax = plt.subplots()
                # ax.plot(time,gr)
                # ax.set_title(k)
                growth_rates[k] = self.est_growth_rate(gr,t=time)
        
        # else:
            

        return growth_rates

    def gen_growth_rate_lib_ts(self):
        """Generates growth rate library from timeseries OD data

        Returns:
            dict: Dict of dose-response curves indexed by replicate
        """
        growth_rates = self.get_growth_rates_from_df()
        replicate_num = 0
        growth_rate_lib = {}

        if self.genotype_dict is not None:
            time = np.array(self.data['Time [s]'])
            for key in self.genotype_dict.keys():
                rate_est = []
                for well in self.genotype_dict[key]:
                    # print(well)
                    ts = np.array(self.data[well])

                    if self.tmax is not None:
                        indx = np.argwhere(time>=self.tmax)
                        indx = indx[0][0]
                        ts = ts[:indx]
                        time_t = time[:indx]

                    d,pcov = self.est_logistic_params(ts,time_t,self.debug,normalize=False)
                    rate_est.append(d['gr']*3600)
                grl_t = {'avg':np.mean(rate_est),
                         'std':np.std(rate_est)}
                growth_rate_lib[key] = grl_t
            return growth_rate_lib

        
        if self.replicate_arrangement == 'rows': # letters represent individual replicates
            # get all the letters in the data keys
            replicates = []
            concentrations = []
            for key in self.data_keys:
                replicates.append(key[0])
                concentrations.append(key[1:])

            replicates = list(set(replicates))
            concentrations = list(set(concentrations))

            replicates.sort()
            
            concentrations = [int(c) for c in concentrations]

            concentrations.sort()
            concentrations = [str(c) for c in concentrations]
            
            for r in replicates:
                # gr_vect = []
                gr_dict = {}
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = r + concentrations[i]
                    gr_dict[str(c)] = growth_rates[key]
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_dict
                replicate_num += 1

            
        else:
            replicates = []
            concentrations = []
            for key in self.data_keys:
                replicates.append(key[1:])
                concentrations.append(key[0])

            replicates = list(set(replicates))
            concentrations = list(set(concentrations))

            concentrations.sort()

            replicates = [int(r) for r in replicates]
            replicates.sort()
            replicates = [str(r) for r in replicates]

            for r in replicates:
                # gr_vect = []
                gr_dict = {}
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = concentrations[i] + r
                    # gr_vect.append(growth_rates[key])
                    gr_dict[str(c)] = growth_rates[key]
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_dict
                replicate_num += 1

        self.n_genotype = replicate_num
        return growth_rate_lib

    def logistic_growth_curve(self,t,r,p0,k):
        """Logistic growth equation

        Args:
            t (float): time
            r (float): growth rate
            p0 (float): starting population size
            k (float): carrying capacity

        Returns:
            float: population size at time t
        """
        p = k/(1+((k-p0)/p0)*np.exp(-r*t))

        return p

    def est_background(self,data=None,background_keys=None):
        if background_keys is None:
            background_keys=self.background_keys
        if data is None:
            data = self.data
        data_dict = self.od_data_to_dict(data)

        bg_est = 0
        for key in background_keys:
            bg_est += data_dict[key]
        bg_est = bg_est/len(background_keys)
        return bg_est

    def set_background(self,data=None,background_keys=None):
        self.background_od = self.est_background(data=data,
                                                 background_keys=background_keys)


    def gen_growth_rate_lib_sm(self):
        """Generate growth rate library from a single measurement experiment

        Returns:
            dict: Dict of dicts. Keys are genotypes. Sub-dicts contain mean growth rate
            ('avg') and standard deviation ('std')
        """
        
        gd = self.parse_exp_layout_file()
        gr_lib = {}
        data_dict = self.od_data_to_dict(self.data)

        for g in gd.keys():
            if g != 'CONTROL':
                keys = gd[g]
                r_t = []
                for k in keys:

                    od = data_dict[k] - self.background_od
                    r = self.OD_rate_eqn(od)

                    if not pd.isna(r):
                        r_t.append(r)

                gr_t = {'avg':np.mean(r_t),
                        'std':np.std(r_t)}
                gr_lib[g] = gr_t

        return gr_lib

    def compute_constants(self,ref_params=None):

        if ref_params is None:
            ref_params = self.ref_params

        OD_max = 0
        OD_0 = 0
        count = 0

        for key in ref_params.keys():
            OD_max += ref_params[key]['OD_max']
            OD_0 += ref_params[key]['OD_0']
            count+=1

        OD_max = OD_max/count
        OD_0 = OD_0/count
        L = (OD_max - OD_0)/OD_0

        constants = {'OD_max':OD_max,
                     'OD_0':OD_0,
                     'L':L}
        return constants

    def OD_rate_eqn(self,OD_obs,t_obs=None,OD_max=None,L=None):
        """Estimates the growth rate using a single OD measurement

        Args:
            t_obs (int): Time of observation (seconds)
            OD_max (float): Estimated max OD
            OD_obs (float): Observed OD
            L (float): Experimental constant ((OD_max - OD_0)/OD_0)

        Returns:
            float: Growth rate (s^-1)
        """
        if t_obs is None:
            t_obs = self.t_obs
        if OD_max is None:
            OD_max = self.OD_max
        if L is None:
            L = self.L

        if OD_obs<=0:
            r = 0
        else:
            r = (-1/t_obs)*np.log((OD_max-OD_obs)/(OD_obs*L))
        
        if r < 0:
            r = 0
        return r
    
    def parse_exp_layout_file(self):
        """Return a dict of genotypes with locations.

        Each entry in the dict contains a list of keys referring to wells in the plate
        where data for the genotype is found.

        Args:
            df (pandas dataframe): Experiment layout dataframe.

        Returns:
            dict: Dict of genotypes where each genotype contains a list of wells
        """
        data_path = self.exp_layout_path
        if '.csv' in data_path:
            df = pd.read_csv(data_path)
        elif '.xlsx' in data_path:
            df = pd.read_excel(data_path)

        # Get the total number of genotypes
        cur_max = 0
        for col in df.columns:
            if col != 'row':
                for d in df[col]:
                    data_type = type(d)
                    if type(d) == int:
                        if d > cur_max:
                            cur_max = d
                    elif d.isnumeric():
                        if int(d) > cur_max:
                            cur_max = int(d)
        
        # For each genotype, generate keys that are associated with that 
        # genotype's location in the plate

        genotype_dict = {}

        for g in np.arange(cur_max+1):

            # generate a list of tuples that are row-col pairs
            if data_type == int:
                indx = df[df == g].stack().index.tolist()
            elif data_type == str:
                indx = df[df == str(g)].stack().index.tolist()

            # transform l into wells (i.e. A1, etc)
            gen_locs = []
            for l in indx:
                # print(l)
                row = l[0]
                col = l[1]

                if not(row == 'row' or col == 'row'):
                    if type(col) == int:
                        key0 = df['row'][row]
                        key1 = col

                    if col.isnumeric(): # rows of the dataframe are letters
                        key0 = df['row'][row]
                        key1 = col

                    else: # rows of the dataframe are numbers
                        key0 = col
                        key1 = df['row'][row]
                    key = key0+str(int(key1))
                    # print(key)
                    gen_locs.append(key)
            genotype_dict[str(g)] = gen_locs
        
        # Get the location of control wells

        indx = df[df=='CONTROL'].stack().index.tolist()
        control_locs = []
        for l in indx:
            row = l[0]
            col = l[1]

            if col.isnumeric(): # rows of the dataframe are letters
                key0 = df['row'][row]
                key1 = col

            else: # rows of the dataframe are numbers
                key0 = col
                key1 = str(int(df['row'][row]))

            key = key0+key1
            control_locs.append(key)
        
        genotype_dict['CONTROL'] = control_locs
        self.background_keys = control_locs

        return genotype_dict

    def parse_od_data_file(self,data_path):
        """Loads the raw OD data files and strips metadata

        OD data files refers to the type of data that is a single OD reading in time (i.e
        not timeseries data).

        Args:
            data_path (str): path to data file

        Returns:
            pandas dataframe: Formatted dataframe with metadata stripped
        """
        if '.csv' in data_path:
            df = pd.read_csv(data_path)
        elif '.xlsx' in data_path:
            df = pd.read_excel(data_path)

        # get the first column as an array
        col_0 = df.columns[0]
        col_0_array = np.array(df[col_0])

        data_start_indx = np.argwhere(col_0_array=='<>')
        data_start_indx = data_start_indx[0][0]

        # the data end index should be the first NaN after the data start index
        col_0_array_bool = [pd.isna(x) for x in col_0_array]
        data_end_indx = np.argwhere(col_0_array_bool[data_start_indx:])
        data_end_indx = data_end_indx[0][0] + data_start_indx - 1

        df_filt = df.loc[data_start_indx:data_end_indx,:]

        # fix the columns
        i = 0
        columns = list(df_filt.iloc[0])
        columns_t = []
        for c in columns:
            if type(c) is not str:
                columns_t.append(str(int(c)))
            else:
                columns_t.append(c)
            i+=1
        
        columns_t[0] = 'Rows'

        df_filt.columns = columns_t

        df_filt = df_filt.drop(df_filt.index[0])
        df_filt = df_filt.reset_index(drop=True)

        return df_filt

    def od_data_to_dict(self,df):
        """Takes an OD data file and returns a dict with each data well as a key

        OD data files refers to the type of data that is a single OD reading in time (i.e
        not timeseries data).

        Args:
            df (pandas dataframe): Parsed OD data file

        Returns:
            dict: Dict of key-value pairs where each key is a well location and each value
            is the OD reading.
        """
        rownum = 0
        d = {}

        for k0 in df['Rows']:
            for k1 in df.columns[1:]:
                if k0.isnumeric():
                    key = k1+k0
                else:
                    key = k0+k1
                d[key] = df[k1][rownum]
            rownum+=1

        return d

    def get_reference_params(self,genotypes=None,keys=None,df=None):
        """Gets the growth rates for one or more wells

        Given a genotype and key or list of genotypes and keys, estimates the growth 
        rate for each genotype using the keys provided.

        Args:
            genotypes (int or list of int): Genotype(s) analyzed
            keys (str or list of str): Plate well keys correspoding to genotypes
            df (pandas dataframe): data

        Returns:
            dict: Dict of results. Keys are genotypes and entries are growth rates.
        """
        if genotypes is None:
            genotypes = self.ref_genotypes
        if keys is None:
            keys = self.ref_keys
        if df is None:
            df = self.ref_data

        time = np.array(df['Time [s]'])
        
        ref_params = {}
        if type(keys) == list:
            g_indx = 0
            for k in keys:

                od = np.array(df[k]) # growth rate time series data
                d = self.est_logistic_params(od,t=time)

                g = genotypes[g_indx]
                ref_params[str(g)] = d

                g_indx +=1
        
        else:
            od = np.array(df[keys])
            d = self.est_logistic_params(od,t=time)
            ref_params[str(genotypes)] = d

        return ref_params

    def est_logistic_params(self,growth_curve,t,debug=False,sigma=None,mode='logistic',
                        normalize=False):
        """Estimates growth rate from OD growth curve

        Args:
            growth_curve (list or numpy array): vector of OD data
            t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

        Returns:
            dict: Dict of logistic growth curve paramters
        """
        # interpolate negative time

        # t_ext = t[0]

        # t = [tt+t[1] for tt in t]

        # t_ext = [t_ext] + t

        # t = t_ext

        # growth_curve = [growth_curve[0]] + growth_curve

        # normalize
        if normalize:
            norm_factor = max(growth_curve)
            growth_curve = [g/norm_factor for g in growth_curve]
        else:
            norm_factor = 1

        # estimate cc

        cc_est = np.mean(growth_curve[-2:])

        gr_est = self.rolling_regression(t,growth_curve)
        
        bounds = ([gr_est-0.2*gr_est,growth_curve[0]-0.1,cc_est-cc_est*0.05,-1000],
                [gr_est+0.2*gr_est,growth_curve[0]+0.1,cc_est+cc_est*0.05,max(t)])

        # p0 = [10**-3,growth_curve[0],cc_est] # starting parameters

        # popt, pcov = sciopt.curve_fit(logistic_growth_curve,
        #                                     t,growth_curve,p0=p0,sigma=sigma,
        #                                     bounds=bounds)
        p0 = [gr_est,growth_curve[0],cc_est,1000]

        popt,pcov = sciopt.curve_fit(self.logistic_growth_with_lag,t,growth_curve,p0=p0,
                                    bounds=bounds)
        
        rate_indx = 0 # index of the optimized data points referring to the growth rate
        p0_indx = 1 # index of the optimized data points referring to initial population size
        carrying_cap_indx = 2 # index of the optmized data points referring to carrying capacity
        lamba_indx = 3

        r = popt[rate_indx]
        p0 = popt[p0_indx]*norm_factor
        cc = popt[carrying_cap_indx]*norm_factor
        l = popt[lamba_indx]

        min_carrying_cap = 0.4

        if r < 0: # if the growth rate is negative
            r = 0
        if cc < p0: # if the carrying capacity is less than the initial population size
            r = 0
        if cc < min_carrying_cap: # if the carrying cap is less the minimum threshold
            r = 0
        if norm_factor < 0.4:
            r = 0
                
        d = {'gr':r,
            'OD_0':p0,
            'OD_max':cc}   

        if debug:
            # if r > 0:
            fig,ax = plt.subplots()
            t_plot = np.array(t)/3600
            ax.scatter(t_plot,growth_curve)

            est = [self.logistic_growth_with_lag(tt,popt[0],popt[1],popt[2],popt[3]) for tt in t]
            
            ax.plot(t_plot,est,color='red')

            p0 = round(popt[1]*10**5)/10**5
            k = round(popt[2]*10**5)/10**5
            r_t = round(popt[0]*3600,2)
            title = 'rate = ' + str(round(3600*r,2)) + ' K = ' + str(k) + ' p0 = ' + str(round(p0,2))

            ax.set_title(title)
            ax.set_xlabel('Time (hr)')
            if normalize:
                ax.set_ylabel('Normalized OD')
            else:
                ax.set_ylabel('OD')
                # ax.set_ylim(0,1)   

        return d,pcov
    
    def logistic_growth_with_lag(self,t,r,p0,k,l):

        p = p0 + k/(1+ np.exp((4*r*(l-t)/k) + 2))

        return p
    
    def rolling_regression(self,xdata,ydata):

        # compute diff

        r = []
        for i in range(len(ydata)-1):
            dy = ydata[i+1] - ydata[i]
            dx = xdata[i+1] - xdata[i]
            r.append(dy/dx)

        cc_est = np.mean(ydata[-2:])

        if cc_est < 0.3:
            return 10**-6
        else:
            return np.max(r)

# Misc helper functions