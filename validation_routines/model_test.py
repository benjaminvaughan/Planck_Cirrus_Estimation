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
hdul = fits.open('/data/butler/SPIRE/hls_clusters/rbs1639_PSW_6_8.2.fits')

# hdul = fits.open('/data/butler/SPIRE/hermes_clusters/rxj1347_PSW_nr_1.fits')
# calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * 300**2

center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
map = hdul[1].data
size = map.shape
ref_head = hdul[1].header

nu = [353e9, 545e9, 857e9]

sim_maps = []
error_maps = []
for i in range(len(nu)):
    I_map = create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=nu[i]) #create_map(filenames, 8, ref_mapsize=ref_size, center=center, nu=nu[i])
    sim_maps.append(I_map)

fig, axs = plt.subplots(1, 3)
im1 = axs[0].imshow(sim_maps[0], origin='lower', vmin=0, vmax=0.32)
fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im2 = axs[1].imshow(sim_maps[1], origin='lower', vmin=0, vmax=1.6)
fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im3 = axs[2].imshow(sim_maps[2], origin='lower', vmin=0, vmax=4)
fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')

plt.tight_layout()
plt.savefig('../Test_Cases/model.png')
plt.clf()


lo_data, pixsize, x_side, y_side, ra, dec= read_in_fits('../Data/COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits', center, ref_head, 6, size)
low_map = np.reshape(lo_data, (x_side, y_side))
interp_low_map = interp_back_to_ref(low_map, ra, dec, ref_head, size)

mi_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits', center, ref_head, 6, size)
mid_map = np.reshape(mi_data, (x_side, y_side))
interp_mid_map = interp_back_to_ref(mid_map, ra, dec, ref_head, size)

hi_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits', center, ref_head, 6, size)
high_map = np.reshape(hi_data, (x_side, y_side))
interp_high_map = interp_back_to_ref(high_map, ra, dec, ref_head, size)

real_maps = [low_map, mid_map, high_map]

for i in range(len(sim_maps)):
    diff_map = real_maps[i] - sim_maps[i]
    print(np.mean(diff_map), i)
    plt.imshow(diff_map)
    plt.colorbar()
    plt.savefig('../Test_Cases/difference_map%2E.png' % nu[i])
    plt.clf()



fig, axs = plt.subplots(1, 3)
im1 = axs[0].imshow(real_maps[0], origin='lower', vmin=0, vmax=0.32)
fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im2 = axs[1].imshow(real_maps[1], origin='lower', vmin=0, vmax=1.6)
fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im3 = axs[2].imshow(real_maps[2], origin='lower', vmin=0, vmax=4)
fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')
plt.tight_layout()
plt.savefig('../Test_Cases/Planck_model.png')
plt.clf()
# plt.clim(0, 15)
# axs[2].set_colorbar().set_label('Flux [$\\frac{MJy}{sr}$]')


lo_data, pixsize, x_side, y_side, ra, dec= read_in_fits('../Data/COM_CompMap_CIB-GNILC-F353_2048_R2.00.fits', center, ref_head, 6, size)
low_map = np.reshape(lo_data, (x_side, y_side))
# interp_low_map = interp_back_to_ref(low_map, ra, dec, ref_head, size)

mi_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F545_2048_R2.00.fits', center, ref_head, 6, size)
mid_map = np.reshape(mi_data, (x_side, y_side))
# interp_mid_map = interp_back_to_ref(mid_map, ra, dec, ref_head, size)

hi_data, pixsize, x_side, y_side, ra, dec = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F857_2048_R2.00.fits', center, ref_head, 6, size)
high_map = np.reshape(hi_data, (x_side, y_side))
# interp_high_map = interp_back_to_ref(high_map, ra, dec, ref_head, size)

cib_maps = [low_map + sim_maps[0], mid_map + sim_maps[1], high_map + sim_maps[2]]

fig, axs = plt.subplots(1, 3)
im1 = axs[0].imshow(cib_maps[0], origin='lower', vmin=0, vmax=0.32)
fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im2 = axs[1].imshow(cib_maps[1], origin='lower', vmin=0, vmax=1.6)
fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im3 = axs[2].imshow(cib_maps[2], origin='lower', vmin=0, vmax=4)
fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')
plt.tight_layout()
plt.savefig('../Test_Cases/CIB_model.png')
plt.clf()
#
sim = []
for i in range(len(sim_maps)):
    sim.append(sim_maps[i] + cib_maps[i])
    plt.imshow(sim[i])
    plt.colorbar()
    plt.savefig('CIBandDUST')
    plt.clf()
    diff_map = real_maps[i] - sim[i]
    print(np.mean(diff_map), i)

# histograms(sim_maps, real_maps, nu)
histograms(sim, real_maps, nu)
