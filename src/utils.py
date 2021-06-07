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
    p = r + os.sep + 'results' + os.sep + filename
    return p