################################################################################
# NAME : pointing_test.py
# DATE STARTED : July 2, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : The purpose of this script is to test the astrometry routines.
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
from make_map import *
from histograms import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world

tau_name = '../Data/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
temp_name = '../Data/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
beta_name = '../Data/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'

filenames = [tau_name, temp_name, beta_name]
hdul = fits.open('../Data/macs2129_PSW_6_8.2.fits')
# hdul = fits.open('/data/butler/SPIRE/hermes_clusters/rxj1347_PSW_nr_1.fits')
calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6) * 18**2

center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
print(center)
map = hdul[1].data
size = map.shape
ref_head = hdul[1].header

x = np.arange(0, size[0])
y = np.arange(0, size[1])
X, Y = np.meshgrid(x, y)

PLW_I_map =  create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=1200e9) * calfac

cib_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F857_2048_R2.00.fits', center, ref_head, 6, size)
CIB_map = np.reshape(cib_data, (x_side, y_side)) * calfac
interped_map = interp_back_to_ref(CIB_map, ra, dec, ref_head, size)

fig, axs = plt.subplots(1, 2)

im1 = axs[0].imshow(map, origin='lower')
fig.colorbar(im1, ax=axs[0])
im2 = axs[1].imshow(PLW_I_map, origin='lower')
fig.colorbar(im2, ax=axs[1])
# im3 = axs[2].imshow(interped_map, origin='lower', vmax=0.032)
# fig.colorbar(im3, ax=axs[2])
plt.tight_layout()
plt.savefig('../Test_Cases/cirrus_test_macs2129.png')
plt.clf()




hdul = fits.open('/data/butler/SPIRE/hls_clusters/rbs1639_PSW_6_8.2.fits')

calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6) * 18**2

# center = [206.9, -11.45]
center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
map = hdul[1].data
size = map.shape
ref_head = hdul[1].header

x = np.arange(0, size[0])
y = np.arange(0, size[1])
X, Y = np.meshgrid(x, y)

PLW_I_map =  create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=1200e9) * calfac

cib_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F857_2048_R2.00.fits', center, ref_head, 6, size)
CIB_map = np.reshape(cib_data, (x_side, y_side)) * calfac
interped_map = interp_back_to_ref(CIB_map, ra, dec, ref_head, size)

fig, axs = plt.subplots(1, 2)

im1 = axs[0].imshow(map, origin='lower')
fig.colorbar(im1, ax=axs[0])
im2 = axs[1].imshow(PLW_I_map, origin='lower')
fig.colorbar(im2, ax=axs[1])
# im3 = axs[2].imshow(interped_map, origin='lower', vmax=0.032)
# fig.colorbar(im3, ax=axs[2])
plt.tight_layout()
plt.savefig('../Test_Cases/cirrus_test_rbs1639.png')
plt.clf()







hdul = fits.open('/data/butler/SPIRE/hls_clusters/macs0451_PSW_6_8.2.fits')

calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6) * 18**2

# center = [206.9, -11.45]
center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
map = hdul[1].data
size = map.shape
ref_head = hdul[1].header

x = np.arange(0, size[0])
y = np.arange(0, size[1])
X, Y = np.meshgrid(x, y)

PLW_I_map =  create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=1200e9) * calfac

cib_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F857_2048_R2.00.fits', center, ref_head, 6, size)
CIB_map = np.reshape(cib_data, (x_side, y_side)) * calfac
interped_map = interp_back_to_ref(CIB_map, ra, dec, ref_head, size)

fig, axs = plt.subplots(1, 2)

im1 = axs[0].imshow(map, origin='lower')
fig.colorbar(im1, ax=axs[0])
im2 = axs[1].imshow(PLW_I_map, origin='lower')
fig.colorbar(im2, ax=axs[1])
# im3 = axs[2].imshow(interped_map, origin='lower', vmax=0.032)
# fig.colorbar(im3, ax=axs[2])
plt.tight_layout()
plt.savefig('../Test_Cases/cirrus_test_macs0451.png')
plt.clf()
