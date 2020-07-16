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
from utilities import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from utilities import make_header
from astropy import units as u

def point(hdul, name):

    tau_name = '../Data/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
    temp_name = '../Data/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
    beta_name = '../Data/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'

    filenames = [tau_name, temp_name, beta_name]

    center = [hdul[1].header['crval1'], hdul[1].header['crval2']]
    map = hdul[1].data
    size = map.shape
    ref_head = hdul[1].header

    x = np.arange(0, size[0])
    y = np.arange(0, size[1])
    X, Y = np.meshgrid(x, y)

    PSW_I_map, ra, dec =  create_map(filenames, ref_head, 6, ref_mapsize=size, center=center, nu=857e9)
    ra  = ra[:,0]
    dec = dec[0, :]
    mid_ra = np.median(ra)
    mid_dec = np.median(dec)
    PSW_header = make_header(6, size, mid_ra, mid_dec)

    hdu = fits.PrimaryHDU(PSW_I_map, PSW_header)
    hdul = fits.HDUList([hdu])
    hdul.writeto('../Test_Cases/new_fits_files/Planck_250_' + name +'.fits', overwrite=True)

    min_dec = np.min(dec)
    max_dec = np.max(dec)
    min_ra  = np.min(ra)
    max_ra  = np.max(ra)

    plt.imshow(PSW_I_map, origin='lower', extent=[min_dec, max_dec, min_ra, max_ra])#, clim=(1.8916812, 8.812404))
    plt.colorbar()
    plt.savefig('../Test_Cases/cirrus_test_' + name + '.png')
    plt.clf()


if __name__ == '__main__':

    DataDir = '../Test_Cases/new_fits_files/'
    hdul = fits.open('../Data/macs2129_PSW_6_8.2.fits')

    point(hdul, 'MACS2129')

    data = hdul[1].data
    hdu = fits.PrimaryHDU(data, hdul[1].header)
    hdul = fits.HDUList([hdu])
    hdul.writeto(DataDir + 'MACS2129.fits', overwrite=True)

    hdul = fits.open('/data/butler/SPIRE/hls_clusters/rbs1639_PSW_6_8.2.fits')

    point(hdul, 'RBS1639')

    data = hdul[1].data
    hdu = fits.PrimaryHDU(data, hdul[1].header)
    hdul = fits.HDUList([hdu])
    hdul.writeto(DataDir + 'RBS1639.fits', overwrite=True)

    hdul = fits.open('/data/butler/SPIRE/hermes_clusters/rxj1347_PSW_nr_1.fits')

    point(hdul, 'RXJ1347')

    data = hdul[1].data
    hdu = fits.PrimaryHDU(data, hdul[1].header)
    hdul = fits.HDUList([hdu])
    hdul.writeto(DataDir + 'RXJ1347.fits', overwrite=True)

    hdul = fits.open('/data/butler/SPIRE/hermes_clusters/ms1054_PSW_nr_1.fits')

    point(hdul, 'MS1054')

    data = hdul[1].data
    hdu = fits.PrimaryHDU(data, hdul[1].header)
    hdul = fits.HDUList([hdu])
    hdul.writeto(DataDir + 'MS1054.fits', overwrite=True)

    hdul = fits.open('/data/butler/SPIRE/hermes_clusters/a0370_PSW_nr_1.fits')

    point(hdul, 'A0370')

    data = hdul[1].data
    hdu = fits.PrimaryHDU(data, hdul[1].header)
    hdul = fits.HDUList([hdu])
    hdul.writeto(DataDir + 'A0370.fits', overwrite=True)
