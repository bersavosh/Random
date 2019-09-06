#! /bin/python
# A small script to estimate uncertainties for low-count Poissonian events
# Based on Gehrels et al. 1986 (https://ui.adsabs.harvard.edu/abs/1986ApJ...303..336G/abstract)
# Written for python 3, but should work in Python 2 as well.
# This script only needs the python Scipy (https://www.scipy.org/) package. 

import argparse
from scipy.stats import poisson
from scipy.optimize import minimize

def lolim(lower, N, cl):
    """
    Function to be used by minimize, to estimate the lower bound for the expected value 
    in a Poisson distribution with N observations.
    
    N: number of events
    
    cl: confidence level to estimate
    
    """
    return abs(poisson.cdf(N - 1, lower) - cl)


def uplim(upper, N, cl):
    """
    Function to be used by minimize, to estimate the upper bound for the expected value 
    in a Poisson distribution with N observations.
    
    N: number of events
    
    cl: confidence level to estimate
    
    """
    return abs(1 - poisson.cdf(N, upper) - cl)


def gehrels(N, cl=0.8413, method='Nelder-Mead' ):
    """
    Estimating Gehrels uncertainties on N observed events. 
    
    N: Number of events
    
    cl: Confidence level to estimate, default value = 0.8413 (1-sigma)
    
    method: minimization method keyword to pass on to scipy.optimize.minimize. 
    In most cases 'Nelder-Mead' (default) seems to work fine.
    
    -- 
    Returns a tuple of lower and upper bounds.
    
    """
    gehrels_lo = minimize(lolim, N, args=(N,cl), method=method)['x'][0]
    gehrels_up = minimize(uplim, N, args=(N,cl), method=method)['x'][0]
    return gehrels_lo, gehrels_up


def parser():
    """
    The parsing function to parse input from the command line.
    """
    parser = argparse.ArgumentParser(description='Script to estimate Gehrels uncertainties on N observed events.')
    parser.add_argument('N', type=int, help='Number of observed events')
    parser.add_argument('--cl', type=float, help='Confidence level to estimate, default value = 0.8413 (1-sigma)')
    args = parser.parse_args()
    return args


##

params = parser()
if params.cl == None:
    params.cl = 0.8413
    results = gehrels(params.N)
else:
    results = gehrels(params.N, params.cl)

print(f'For {params.N} events, {params.cl} confidence level boundaries are:\n{results[0]:.2f} -- {results[1]:.2f}')
