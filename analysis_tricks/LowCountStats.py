#! /bin/python
# A small script to estimate uncertainties for low-count Poissonian events
# Based on 
# 1- Gehrels et al. 1986 (https://ui.adsabs.harvard.edu/abs/1986ApJ...303..336G/abstract)
# 2- Kraft et al. 1991 (https://ui.adsabs.harvard.edu/abs/1991ApJ...374..344K/abstract)
# Works with Python 2.6 and higher.
# This script only needs the python Numpy (https://numpy.org/) and Scipy (https://www.scipy.org/) packages. 

from __future__ import print_function
import argparse
import numpy as np
from math import factorial
from scipy.stats import poisson
from scipy.optimize import minimize
from scipy.integrate import quad

### Functions for Gehrels method

def gehrels_lolim(lower, N, cl):
    """
    Function to be used by minimize, to estimate the lower bound for the expected value 
    in a Poisson distribution with N observations.
    
    N: number of events
    
    cl: confidence level to estimate
    
    """
    return abs(poisson.cdf(N - 1, lower) - cl)


def gehrels_uplim(upper, N, cl):
    """
    Function to be used by minimize, to estimate the upper bound for the expected value 
    in a Poisson distribution with N observations.
    
    N: number of events
    
    cl: confidence level to estimate
    
    """
    return abs(1 - poisson.cdf(N, upper) - cl)


def gehrels(N, cl=0.8413, method='Nelder-Mead'):
    """
    Estimating Gehrels uncertainties on N observed events. 
    
    N: Number of events
    
    cl: Confidence level to estimate, default value = 0.8413 (1-sigma)
    
    method: minimization method keyword to pass on to scipy.optimize.minimize. 
    In most cases 'Nelder-Mead' (default) seems to work fine.
    
    -- 
    Returns a tuple of lower and upper bounds.
    
    """
    gehrels_lo = minimize(gehrels_lolim, N, args=(N,cl), method=method)['x'][0]
    gehrels_up = minimize(gehrels_uplim, N, args=(N,cl), method=method)['x'][0]
    return gehrels_lo, gehrels_up

### Functions for Kraft method

def kraft_posterior(S, N, B):
    cn = 0
    for n in np.arange(0,N+1):
        cn += (np.exp(-B) * B**n)/factorial(n)
    C = cn**-1
    return C * (np.exp(- S - B) * (S + B)**N) / factorial(N)


def kraft_posterior_root(S, N, B, Y):
    return abs(kraft_posterior(S, N, B) - Y)


def kraft_posterior_prob(Y, N, B):
    lolim = minimize(kraft_posterior_root, N-B-1, args=(N,B,Y),method='Nelder-Mead')['x'][0]
    if lolim < 0.0:
        lolim = 0.0
    uplim = minimize(kraft_posterior_root, N+B+1, args=(N,B,Y),method='Nelder-Mead')['x'][0]
    if uplim < 0.0:
        uplim = 0.0
    return lolim, uplim, quad(kraft_posterior, lolim, uplim, args=(N,B))[0]


def kraft_post_prob_root(Y, N, B, cl):
    return abs(kraft_posterior_prob(Y, N, B)[2] - cl)


def kraft(N, B, cl):
    starting_point = 0.5*kraft_posterior(N-B, N, B)
    return kraft_posterior_prob(minimize(kraft_post_prob_root, starting_point, (N, B, cl), method='Nelder-Mead')['x'], N, B)


###

def parser():
    """
    The parsing function to parse input from the command line.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Script to estimate uncertainties in low-count statistics based on \
Gehrels et al. 1986, ApJ, 303, 336 and Kraft et al. ApJ 374, 344. \
If no information about the background is provided, the estimated \
limits will be based on Poisson statistics as prescribed by Gehrels. \
Otherwise, they will be limits on the background-subtracted number \
of events, following the prescription by Kraft.\n\n\
\
Examples:\n\n\
\
Estimating Gehrels 1-sigma confidence limits for 6 counts:\n\
> python LowCountStats.py 6\n\n\
\
Estimating Gehrels 0.95 confidence limits for 5 counts:\n\
> python LowCountStats.py 5 --cl 0.95\n\n\
\
Estimating Kraft 1-sigma confidence limits for 7 counts and 4.2 bkg counts\n\
> python LowCountStats.py 7 --bkg 4.2\n\n\
\
Estimating Kraft 0.95 confidence limits for 11 counts and 6 bkg counts\n\
> python LowCountStats.py 11 --bkg 6 --cl 0.95.\n\n\
--------------') 
    parser.add_argument('N', type=int, help='Number of observed events')
    parser.add_argument('--cl', type=float, help='Confidence level to estimate, default value = 0.8413 (1-sigma)')
    parser.add_argument('--bkg', type=float, help='Scaled background counts. If not provided, limits are estimated \
based on Poisson statistics as in Gehrels+86. If provided, the \
background-subtracted limits will be estimated based on the \
Bayesian posterior provided by Kraft+91.\n\n')
    args = parser.parse_args()
    return args


##

params = parser()
if params.cl == None:
    params.cl = 0.8413

if params.bkg == None:
    results = gehrels(params.N, params.cl)
    #print(f'For {params.N} events, Gehrels {params.cl} confidence level boundaries are:\n{results[0]:.2f} -- {results[1]:.2f}')
    print('For '+str(params.N)+' events, Gehrels '+str(params.cl)+' confidence level boundaries are:\n')
    print(round(results[0],2),'--',round(results[1],2),'\n')
    #print('\nReference: Gehrels et al. 1986, ApJ, 303, 336')
    
else:
    results = kraft(params.N, params.bkg, params.cl)
    #print(f'For {params.N} events and a background of {params.bkg}, Kraft {params.cl} confidence level boundaries are:\n{results[0]:.2f} -- {results[1]:.2f}')
    print('For '+str(params.N)+' events and a background of '+str(params.bkg)+', Kraft '+str(params.cl)+' confidence level boundaries are:\n')
    print(round(results[0],2),'--',round(results[1],2),'\n')
    #print('\nReference: Gehrels et al. 1986, ApJ, 303, 336; Kraft et al. 1991, ApJ, 374, 344')
