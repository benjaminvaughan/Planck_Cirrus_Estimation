################################################################################
# NAME : model_test.py
# DATE STARTED : July 5, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : The purpose of this script is to test that the map we recreate using
# the best fit params from Planck agree with their model.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
from astropy.io import fits
sys.path.append('../')
from utilities import *
from make_map import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world

#this is a reference image
hdul = fits.open('/data/butler/SPIRE/hermes_clusters/ms1054_PSW_nr_1.fits')

name = 'ms1054'
ref_head = hdul[1].header
pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])

PSW_I_map, ra, dec =  create_map(ref_head, nu=1200e9)
# PSW_I_map = np.multiply(PSW_I_map, calfac)
ra  = ra[:,0]
dec = dec[0, :]
mid_ra = np.median(ra)
mid_dec = np.median(dec)
PSW_header = make_header(pixsize, PSW_I_map.shape, mid_ra, mid_dec)
hdu = fits.PrimaryHDU(PSW_I_map, PSW_header)
hdul = fits.HDUList([hdu])
hdul.writeto('../Test_Cases/fits_files/Planck_250_' + name +'.fits', overwrite=True)

min_dec = np.min(dec)
max_dec = np.max(dec)
min_ra  = np.min(ra)
max_ra  = np.max(ra)

plt.imshow(PSW_I_map, origin='lower', extent=[min_dec, max_dec, min_ra, max_ra])#, clim=(1.8916812, 8.812404))
plt.colorbar()
plt.savefig('../Test_Cases/' + name + '.png')
plt.clf()
