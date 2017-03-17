from __future__ import with_statement

def readbins(filename):
    """
        read a set of bins in Carlo Contaldi format:
         binmin binmax bintype
         with bintype = 1, 2, 3 for TT, TE, EE
         
         output is a list of lists 
         [ 
            [(TTstart1, TTend1), (TTstart2, TTend2), ...]
            [(TEstart1, TEend1), (TEstart2, TEend2), ...]
            [(EEstart1, EEend1), (EEstart2, EEend2), ...]
         ]
    """
    
    lmin = []
    lmax = []
    
    binlist = [ [], [], [] ]
    with file(filename) as f:
        for line in f:
            lmin, lmax, bintype = [int(i) for i in line.split()]
            binlist[bintype-1].append( (lmin, lmax) )
            
    return binlist