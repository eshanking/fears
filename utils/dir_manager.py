from pathlib import Path
import os

def get_project_root() -> Path:
    return Path(__file__).parent.parent

def make_datapath_absolute(filename):
    r = str(get_project_root())    
    p = r + os.sep + 'data' + os.sep + filename
    return p

def make_resultspath_absolute(filename):
    r = str(get_project_root())    
    res_dir = r + os.sep + 'results'
    make_directory(res_dir)
    p = r + os.sep + 'results' + os.sep + filename
    
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