################################################################################
# NAME : make_map.py
# DATE STARTED : July 2, 2020
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
################################################################################import sys
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
hdul = fits.open('/data/butler/SPIRE/hermes_clusters/rxj1347_PSW_nr_1.fits')
calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6) * 18**2

print(hdul[1].header)

hdul.info()
# center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
center = [206.9, -11.45]
map = hdul[1].data
size = map.shape
ref_head = hdul[1].header

x = np.arange(0, size[0])
y = np.arange(0, size[1])
X, Y = np.meshgrid(x, y)

w = world(ref_head)
coord = pixel_to_skycoord(X, Y, wcs=w, origin=0)
RA = np.asarray(coord.ra.to_string(decimal=True), dtype='float')
DEC = np.asarray(coord.dec.to_string(decimal=True), dtype='float')

PLW_I_map =  create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=1200e9) * calfac

fig, axs = plt.subplots(1, 2)

im1 = axs[0].imshow(map, origin='lower')
fig.colorbar(im1, ax=axs[0])
im2 = axs[1].imshow(PLW_I_map, origin='lower')
fig.colorbar(im2, ax=axs[1])
plt.savefig('../Test_Cases/cirrus_test_macs2129.png')
plt.clf()


exit()
ref_size = [250, 250]
nu = [353e9, 545e9, 857e9]


sim_maps = []
error_maps = []
for i in range(len(nu)):
    I_map =  create_map(filenames, 8, ref_mapsize=ref_size, center=center, nu=nu[i])
    sim_maps.append(I_map)
    # error_maps.append(E_map)
# convolved_I_map = convolve(I_map, beam, normalize_kernel=False)

vmin = 0
# vmax = 15

lo_data, pixsize= read_in_fits('../Data/COM_CompMap_CIB-GNILC-F353_2048_R2.00.fits', center, 8, ref_mapsize=ref_size)
side = int(sqrt(len(lo_data)))
low_map = np.reshape(lo_data, (side, side))

plt.imshow(low_map)
plt.savefig('../Data/TEST.png')
plt.clf()
exit()
mi_data, pixsize = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F545_2048_R2.00.fits', center,  8, ref_mapsize=ref_size)
side = int(sqrt(len(mi_data)))
mid_map = np.reshape(mi_data, (side, side))

# hi_data, pixsize = read_in_fits('../Data/COM_CompMap_CIB-GNILC-F857_2048_R2.00.fits', center, ref_pixsize=8, ref_mapsize=ref_size)
# side = int(sqrt(len(hi_data)))
# high_map = np.reshape(hi_data, (side,side))


fig, axs = plt.subplots(1, 3)
im1 = axs[0].imshow(sim_maps[0], origin='lower', vmin=0, vmax=3.0)
fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im2 = axs[1].imshow(sim_maps[1], origin='lower', vmin=0, vmax=10.5)
fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im3 = axs[2].imshow(sim_maps[2], origin='lower', vmin=0, vmax=30)
fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')

plt.tight_layout()
plt.savefig('../Test_Cases/model.png')
plt.clf()


# fig, axs = plt.subplots(1, 3)
# im1 = axs[0].imshow(error_maps[0], origin='lower', vmax= (i + 1), vmin=0)
# fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
# im2 = axs[1].imshow(error_maps[1], origin='lower', vmax= (i + 1), vmin=0)
# fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
# im3 = axs[2].imshow(error_maps[2], origin='lower', vmax= (i + 1), vmin=0)
# fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')
#
# plt.tight_layout()
# plt.savefig('error_maps.png')
# plt.clf()


lo_data, pixsize= read_in_fits('../Data/COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits', center, 8, ref_mapsize=ref_size)
side = int(sqrt(len(lo_data)))
low_map = np.reshape(lo_data, (side, side))

mi_data, pixsize = read_in_fits('../Data/COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits', center, 8, ref_mapsize=ref_size)
side = int(sqrt(len(mi_data)))
mid_map = np.reshape(mi_data, (side, side))

hi_data, pixsize = read_in_fits('../Data/COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits', center,  8, ref_mapsize=ref_size)
side = int(sqrt(len(hi_data)))
high_map = np.reshape(hi_data, (side,side))

real_maps = [low_map, mid_map, high_map]

for i in range(len(sim_maps)):
    diff_map = sim_maps[i] - real_maps[i]
    print(np.mean(diff_map), i)
    plt.imshow(diff_map, vmin=-1 *(i + 1), vmax=0)
    plt.colorbar()
    plt.savefig('../Test_Cases/difference_map%2E.png' % nu[i])
    plt.clf()



fig, axs = plt.subplots(1, 3)
im1 = axs[0].imshow(low_map, origin='lower', vmin=0, vmax=3.0)
fig.colorbar(im1, ax=axs[0])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im2 = axs[1].imshow(mid_map, origin='lower', vmin=0, vmax=10.5)
fig.colorbar(im2, ax=axs[1])#.set_label('Flux [$\\frac{MJy}{sr}$]')
im3 = axs[2].imshow(high_map, origin='lower', vmin=0, vmax=30)
fig.colorbar(im3, ax=axs[2])#.set_label('Flux [$\\frac{MJy}{sr}$]')
plt.tight_layout()
plt.savefig('../Test_Cases/Planck_model.png')
plt.clf()
# plt.clim(0, 15)
# axs[2].set_colorbar().set_label('Flux [$\\frac{MJy}{sr}$]')
histograms(sim_maps, real_maps, nu)
