#!/usr/bin/python
# -*- coding: utf-8 -*-
# File: gkfit.py
# Created: 2016-03-09 by gks 
"""
Description: Generalized normal distribution. Used in Top hat fitting for diffusers.

Notebook: http://localhost:8888/notebooks/Top%20-%20Hat%20fitting.ipynb#
"""

from scipy.special import gamma as Gamma
import numpy as np
from scipy.interpolate import UnivariateSpline
import scipy

def helpme():
    print("""
    # Initial parameters
    params = Parameters()
    params.add('C_off', value=np.min(y),vary=False)
    params.add('C_sc', value=200.)
    params.add('mu', value=200.)
    params.add('alpha', value=10.0)
    params.add('beta', value=20.0,min=0)

    # Fitting
    out1 = minimize(gkfit.gen_normal_dist_res_func, params, args=(x, y))
    mymodel1 = y + out1.residual
    lmfit.report_fit(out1)
    """
    )
        
def gen_normal_dist(x,C_off=0.,C_sc=1.,mu=0.,alpha=1.0,beta=2.0):
    """
    Calculates a generalized normal distribution function. If beta = 2, then normal. 
    If you increase beta it gets more top-hatty.
    Initialized for a regular normal distribution.

    INPUTS:
        x      np. array
        C_off  offset constant
        C_sc   scaling constant
        mu     mean
        beta   shape
        alpha  scale (or ~sigma)

    RETURNS:
        Generalized normal

    EXAMPLE:
        x = np.linspace(-3,3,100)
        y = gen_normal_dist(x,mu=0,alpha=1,beta=5.)
        plt.plot(x,y)

    NOTES:
        See PDF here:
        https://en.wikipedia.org/wiki/Generalized_normal_distribution
    """
    return C_sc * ( (beta)/(2.*alpha*Gamma(1./beta)) ) * np.exp(- ( np.abs(x-mu)/alpha )**beta ) + C_off

# Define fitting function
def gen_normal_dist_res_func(params, x, y, eps=None):
    """
    The residual function for gen_normal_dist() to maximize.

    INPUTS:
        params is a Parameter() object from lmfit
        x is a numpy array of x values
        y is a numpy array of y values
        eps is a numpy array with y value errors
    
    OUTPUT:
        residuals to maximize

    EXAMPLE:

    NOTES:
    """
    C_off = params['C_off'].value
    C_sc  = params['C_sc'].value
    mu    = params['mu'].value
    alpha = params['alpha'].value
    beta  = params['beta'].value

    model = gen_normal_dist(x,C_off,C_sc,mu,alpha,beta)

    if eps==None:
        return (model-y)
    else:
        return (model-y)/eps

def calc_fwhm(x,y,return_roots=False):
    """
    Calculates the FWHM for a given dataset, by finding the roots of splines.
    
    INPUTS:
        x - x input array
        y - y input array
    
    OUTPUT:
        FWHM The Full Width Half Max of the data
    
    OPTIONAL:
        return_roots - Boolean. If true, return the roots.
        
    EXAMPLE:
        
        
    """
    spline = UnivariateSpline(x, y-np.max(y)/2.0, s=0)
    roots = spline.roots()
    fwhm = roots[1]-roots[0]
    if return_roots == True:
        return fwhm, roots
    else:
        return fwhm

def calc_fwzm(x,y,z=20,return_roots=False):
    """
    Calculates the Full with at the Z-th maximum for a given dataset, by finding the roots of splines.
    
    INPUTS:
        x - x input array
        y - y input array
    
    OUTPUT:
        FWZM The Full Width Half Max of the data
    
    OPTIONAL:
        return_roots - Boolean. If true, return the roots.
        
    EXAMPLE:
        
        
    """
    try:
        spline = UnivariateSpline(x, y-np.max(y)/z, s=0)
        roots = spline.roots()
        fwzm = roots[1]-roots[0]
        if return_roots == True:
            return fwzm, roots
        else:
            return fwzm
    except Exception as e:
        print("Error with finding roots! Returning 0")
        if return_roots == True:
            return 0., 0
        else:
            return 0.

def calc_hwzm(x,y,z=20):
    """
    Calculates the HWHM at the Z-th maximum for a given dataset, by finding the roots of splines.
    
    INPUTS:
        x - x input array
        y - y input array
    
    OUTPUT:
        HWZM The Half Width at Z-th Max of the data
    
    EXAMPLE:
    """
    spline = UnivariateSpline(x, y-np.max(y)/z, s=0)
    roots = spline.roots()
    return roots

def calc_hwhm(x,y):
    """
    Calculates the HWHM for a given dataset, by finding the roots of splines.
    
    INPUTS:
        x - x input array
        y - y input array
    
    OUTPUT:
        HWZM The Half Width at Z-th Max of the data
    
    EXAMPLE:
    """
    roots = calc_hwzm(x,y,z=2)
    return roots

def line(x,a=1.,b=1.):
    return a*x + b

def line_res_func(params,x,y,eps=None):
    """
    The residual function for line() to maximize.

    INPUTS:
        params is a Parameter() object from lmfit
        x is a numpy array of x values
        y is a numpy array of y values
        eps is a numpy array with y value errors
    
    OUTPUT:
        residuals to maximize

    EXAMPLE:

    NOTES:
    """
    a = params['a'].value
    b = params['b'].value
    
    model = line(x,a,b)
    
    if eps==None:
        return (model-y)
    else:
        return (model-y)/eps


def sinusoid(t,A,P,phi,offset):
    """
    A simple sinusoid model
    """
    return A*np.sin(2.*np.pi*t/P + phi)+offset

def res_func_photometry(params,detrend_array,rel_flux_array):
    """
    A function to detrend photometry

    INPUT:
        params         - detrending coeffs
        detrend_array  - around 1. (normalized)
        rel_flux_array - around 1. (normalized)

    OUTPUT:
        residual detrended array
    
    EXAMPLE:
        residual = res_func_photometry([-0.2,-0.0366],detrend_array,rel_flux_array)
    """
    model = np.zeros(len(rel_flux_array))
    for i,mypar in enumerate(params):
        model += mypar*(detrend_array[i]-1.)
    return (rel_flux_array-1.) - model

def detrend_phot(init_param,detrend_array,rel_flux_array):
    """
    A function to detrend photometry
    
    INPUT:
        init_param = [ a list of initial guess parameters, same length as detrend_array ]
        detrend_array = [ a list of arrays of normalized arrays to detrend with  ]
        rel_flux_array - normalized array of flux values to detrend
        
    RETURNS:
        return residual, p_opt[0], std, chi2
        
    """
    p_opt = scipy.optimize.leastsq(res_func_photometry,init_param,args=(detrend_array,rel_flux_array))
    
    residual = res_func_photometry(p_opt[0],detrend_array,rel_flux_array)
    chi2 = sum(residual*residual)
    std  = np.std(residual)
    
    print("Fit finished.")
    print("Best parameters are",p_opt[0])
    return residual, p_opt[0], std, chi2
