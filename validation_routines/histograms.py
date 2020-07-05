################################################################################
# NAME : validation routines
# DATE STARTED : July 2, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This script contains functions that perform test cases.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import matplotlib.pyplot as plt

def histograms(sim_data, real_data, nu_list):

    for i in range(len(nu_list)):
        sim = sim_data[i].flatten() - np.mean(sim_data[i])
        real = real_data[i].flatten() - np.mean(real_data[i])
        plt.hist(sim, int(50), label='sim', histtype='step')
        plt.hist(real, int(50), label='real' , histtype='step')
        plt.xlabel('Flux [$\\frac{MJy}{Sr}$]')
        plt.ylabel('Nsources')
        plt.legend()
        plt.savefig('../Test_Cases/validation histograms %s.png' % (nu_list[i]))
        plt.clf()
