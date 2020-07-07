################################################################################
# NAME : make_map.py
# DATE STARTED : June 25, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This set of code contains the mathematical functions for calculating
# cirrus emission using a modified blackbody from https://arxiv.org/pdf/1312.1300.pdf
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
from math import *

def blackbody_func(T, nu):
    k       = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    c       = 2.99792458e8            #m/s
    const =  2 * h * nu**3 / c**2
    exponent = np.divide(h * nu / k, T)

    denom = np.exp(exponent) - 1
    bb = np.divide(const, denom)
    return bb

def calc_intensity(beta, tau,T, nu):
    nu_const = nu / 353e9 #GHz
    print(nu_const)
    power = np.asarray([nu_const**bi for bi in beta])
    bb = blackbody_func(T, nu)
    I = np.multiply(tau, bb)
    I_nu = np.multiply(I, power)
    return I_nu
