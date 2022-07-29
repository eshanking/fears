from pathlib import Path
import os
import pandas as pd
import numpy as np
import pickle

def get_project_root() -> Path:
    return Path(__file__).parent.parent

def make_datapath_absolute(filename):
    r = str(get_project_root())    
    p = r + os.sep + 'data' + os.sep + filename
    return p

def make_resultspath_absolute(filename):
    r = str(get_project_root())
    
    # check if filename is already absolute
    if r not in filename:
        
        res_dir = r + os.sep + 'results'
        make_directory(res_dir)
        p = r + os.sep + 'results' + os.sep + filename
    else:
        p = filename
    
    return p

def make_figurepath_absolute(filename):
    r = str(get_project_root())  
    fig_dir = r + os.sep + 'figures'
    make_directory(fig_dir)
    p = r + os.sep + 'figures' + os.sep + filename
    return p

def make_directory(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)
        
    # Load data
    # also use to load ic50 and drugless growth rate (anything from a csv)
def load_fitness(data_path):
    fitness = pd.read_csv(data_path)
    cols = list(fitness.columns)
    fit_array = np.array(cols)
    fit_array = fit_array.astype('float')
    return fit_array

def load_growth_rate_data(data_path):
    
    data = pd.read_csv(data_path)
    
    return data

def load_experiment(exp_path):
    e = pickle.load(exp_path)
    return e
