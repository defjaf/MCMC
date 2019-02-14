

import string

from . import ClData_CosmoMC
from . import ClData_WMAP
from . import ClData_bcp

def ClData(file):
    """
    factory function which returns an instance of the appropriate ClData
    class depending on the file name
    """
    
    if file == 'MAP' or file == 'WMAP':
        return ClData_WMAP.ClData_WMAP()
    elif (string.split(file, '.')[-1]=='newdat'):
        return ClData_bcp.ClData_bcp(file)
    else:
        return ClData_CosmoMC.ClData_CosmoMC(file)

def getClData(listfile, verbose=True, no_pol=False):
    """
    get a list of (CosmoMC format) ClData objects as listed in the file
    blank lines and those beginning with '#' are ignored
    """
    
    data = []
    for dataset in file(listfile):
        dataset=dataset.strip()
        if len(dataset) and dataset[0] != "#":
            if verbose: print("Getting %s" % (dataset))
            set = ClData(dataset)
            if no_pol:
                set.has_pol_really = set.has_pol
                set.has_pol=False   # explicitly ignore polarization
            if verbose: print(" got %s" % (set.name))
            data.append(set)
    return data
