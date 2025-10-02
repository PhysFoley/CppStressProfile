# Functions for performing simple blocking analysis of correlated time series
# data. Based on the method presented in:
#     H. Flyvbjerg and H.G. Petersen, J. Chem. Phys. _91_, 461 (1989)

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

# Theoretical blocking curve for exponentially correlated data series
# with a single decay time, t
def theory(k, t):
    c = np.exp(-1./t)
    b = np.power(2.,k)
    num = c*(1.-np.power(c,b)-b+(c*b))
    den = b*np.power(c-1.,2)
    
    return np.sqrt(1. - (2.*num/den) )

# Block data in given python list
#
# Arguments:
#   data : python list of floats or numpy array
#   verbose : controls how much text the function outputs
# Returns:
#   tuple (mu, c, sig, dsig)
#       mu : mean of data
#       c : return value of curve_fit (2-tuple of optimal parameters and param covariance matrix), c=False if fit failed
#       sig : a list with the successive estimated error on the mean at each blocking step
#       dsig : a list with estimated errors on sig
def block(data, verbose=True):
    data = np.array(data)
    sig  = []
    dsig = []
    
    #calculate estimator of mean
    mu = data.mean()
    
    if verbose: print("\nComputing error estimates on continually blocked data...")
    
    while len(data) >= 2:
        sig.append( data.std(ddof=1)/math.sqrt(float(len(data))) )
        dsig.append( sig[-1] / math.sqrt( 2.0 * (float(len(data)) - 1.0) ) )
        
        if verbose: print("datapoints: {:<5d}mean: {:<17.10f} error: {:<17.10f}".format(len(data), m, sig[-1]) )
        
        # blocking averaging:
        data = np.array([0.5*(x[0]+x[1]) for x in zip(data[0::2],data[1::2])])
    
    #normalize for fitting
    s = sig/sig[0]
    ds = dsig/sig[0]
    
    initguess = [20.]
    try:
        c = scipy.optimize.curve_fit(theory, range(len(s)), s, sigma=ds, absolute_sigma=True, p0=initguess, method="trf", bounds=( [0.],[np.inf]) )
    except RuntimeError:
        c = False
        print("Blocking curve fit failed")
    
    return (mu, c, sig, dsig)
