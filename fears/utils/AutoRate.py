import pandas as pd
import os
import scipy.optimize as sciopt
import scipy.interpolate as sciinter
import matplotlib.pyplot as plt
import numpy as np

class Experiment():
    """Experiment class for a given plate reader experiment
    """
    def __init__(self,
                folder_path,
                moat=False,
                replicate_arrangement='rows',
                drug_conc = None,
                units = 'ug/mL',
                data_cols = None,
                debug=False):
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
            p = Plate(pdp,self.drug_conc,debug=self.debug,moat=self.moat,replicate_arrangement=self.replicate_arrangement,data_cols=dc)
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
        offset = 0
        for p in self.plates:

            keys = [int(k) for k in p.growth_rate_lib.keys()]
            
            for k in keys:
                rep_num =  k + offset
                gl[str(rep_num)] = p.growth_rate_lib[str(k)]
            
            offset += max(keys) + 1

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

        replicates = [int(k) for k in self.growth_rate_lib.keys()]
        
        for r in replicates:

            popt = self.fit_hill_curve(dc,gl[str(r)])
            
            # inidices of optimized parameter vector
            ic50_indx = 0
            g_drugless_indx = 1
            hill_coeff_indx = 2

            ic50 = popt[ic50_indx]
            g_drugless = popt[g_drugless_indx]
            hill_coeff = popt[hill_coeff_indx]

            d_t = {'ic50':ic50,
                'g_drugless':g_drugless,
                'hill_coeff':hill_coeff}

            sl[str(r)] = d_t

        return sl

    def fit_hill_curve(self,xdata,ydata):
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

        
        p0 = [0,ydata[0],-0.08]

        if ydata[0] == 0:
            g_drugless_bound = [0,1]
        else:
            # want the estimated drugless growth rate to be very close to the value given in ydata
            g_drugless_bound = [ydata[0]-0.0001*ydata[0],ydata[0]+0.0001*ydata[0]]

        bounds = ([-5,g_drugless_bound[0],-1],[4,g_drugless_bound[1],-0.001]) # these aren't magic numbers these are just starting parameters that happen to work

        popt, pcov = sciopt.curve_fit(self.logistic_pharm_curve_vectorized,
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


class Plate():
    """96-well plate object
    """
    def __init__(self,
                 data_path,
                 drug_conc,
                 replicate_arrangement='rows',
                 moat=False,
                 debug=False,
                 data_cols=None):
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
        self.data = self.parse_data_file(data_path)
        self.data_cols = data_cols

        if drug_conc is None:
            self.drug_conc = [0,0.003,0.0179,0.1072,0.643,3.858,23.1481,138.8889,833.3333,5000]
        else:
            self.drug_conc = drug_conc

        self.replicate_arrangement = replicate_arrangement
        self.debug = debug

    def execute(self):

        self.background_keys = self.get_background_keys()
        self.data_keys = self.get_data_keys()
        self.growth_rate_lib = self.gen_growth_rate_lib()

    def parse_data_file(self,p):
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
        
        # raw data starts after cycle nr.
        if any(df.keys() == 'Cycle Nr.'):
            data_start_indx = np.argwhere(list(time_array) == 'Cycle Nr.')
        elif any(df.keys() == 'Time [s]'):
            data_start_indx = np.argwhere(list(time_array) == 'Time [s]')
        else:
            raise Exception('Unknown file format. Expected either Cycle Nr. or Time [s] as column headings.')

        #sometimes the data gets passed in very unraw
        if len(data_start_indx) == 0:
            return df
        
        data_start_indx = data_start_indx[0][0] # get scalar from numpy array

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
    
    def est_growth_rate(self,growth_curve,t=None,carrying_cap=4):
        """Estimates growth rate from OD growth curve

        Args:
            growth_curve (list or numpy array): vector of OD data
            t (list or numpy array, optional): Time vector. If None, algorithm assumes each time step is 1 s. Defaults to None.

        Returns:
            float: Growth rate in units of 1/s
        """
        
        if t is None:
            t = np.arange(len(growth_curve))

        growth_curve = growth_curve/carrying_cap

        p0 = [10**-6,0.05,1] # starting parameters

        popt, pcov = sciopt.curve_fit(self.logistic_growth_curve,
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

        data_keys = self.get_data_keys()
        time = df['Time [s]']

        for k in data_keys:
            gr = np.array(df[k]) # growth rate time series data
            time = np.array(time) # time vector
            growth_rates[k] = self.est_growth_rate(gr,t=time)

        return growth_rates

    def gen_growth_rate_lib(self):
        """Generates growth rate library from OD data

        Returns:
            dict: Dict of dose-response curves indexed by replicate
        """
        growth_rates = self.get_growth_rates_from_df()
        replicate_num = 0
        growth_rate_lib = {}

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
            concentrations.sort()
            
            for r in replicates:
                gr_vect = []
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = r + concentrations[i]
                    gr_vect.append(growth_rates[key])
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_vect
                replicate_num += 1

        else:
            replicates = []
            concentrations = []
            for key in self.data_keys:
                replicates.append(key[1:])
                concentrations.append(key[0])

            replicates = list(set(replicates))
            concentrations = list(set(concentrations))

            replicates.sort()
            concentrations.sort()

            for r in replicates:
                gr_vect = []
                i = 0 # concentration index
                for c in self.drug_conc:
                    key = concentrations[i] + r
                    gr_vect.append(growth_rates[key])
                    i += 1

                growth_rate_lib[str(replicate_num)] = gr_vect
                replicate_num += 1

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