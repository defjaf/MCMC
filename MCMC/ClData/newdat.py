

def write_newdat():
    """ write C_l data in *.newdat format 
    
    Read a dataset in the following format from Carlo Contaldi
    (courtesy CosmoMC)
    !name
    !n_bands_TT n_EE, n_BB, n_EB, n_TE, n_TB
    !has_calib_uncertain calib calib_err
    !has_beam_uncertain beam beam_err
    !ilike (0: Gaussian, 1: all x-factor, 2: specified have x-factor)
    !loop over {
    ! band-types
    ! band info: num obs + - x_factor l_min l_max use_x
    ! correlation-matrix (ignored)
    ! }  
    ! covariance matrix
    
    """
    
    pass
    
def readraw(prefix):
    """ read the datafiles written by test_binnedCl:
    prefix.bp      ell, mean, err1, err2
    prefix.bpte
    prefix.bpee

    prefix.covar    [savetxt format]
    
    prefix.corr
    
    prefix(nn)      l, Wl[0], Wl[1], Wl[2]
    
    need to write/read xb files
    
    
    """
    
    pass
    