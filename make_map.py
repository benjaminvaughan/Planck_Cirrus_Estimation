################################################################################
# NAME : make_map.py
# DATE STARTED : June 25, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE :
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from astropy.io import fits
import numpy as np
# import healpy.pixelfunc as hp
from math import *
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve_fft as convolve
from astropy.convolution import Gaussian2DKernel

def create_map(filenames, ref_pixsize, ref_mapsize, center):
    param_values = [] #0 = Tau, 1 = Temp, 2 = Emissivity
    for f in filenames:
        data, pixsize = read_in_fits(f, center, ref_pixsize=8, ref_mapsize=ref_size)
        side = int(sqrt(len(data)))
        param_values.append(data)
    nu = 353*1e9 #beta, tau,T, nu
    I = calc_intensity(param_values[2], param_values[0], param_values[1], nu)

    beam, MJy_conv = get_fwhm(center, ref_pixsize=pixsize, ref_mapsize=ref_size)

    I_map = np.reshape(I, (side, side)) * MJy_conv
    return I_map

def get_fwhm(center, ref_pixsize=8, ref_mapsize=260):
    hdul_fwhm = fits.open('COM_CompMap_Dust-GNILC-Beam-FWHM_0128_R2.00.fits')

    header = hdul_fwhm[1].header
    nside = header['NSIDE']
    order = header['ORDERING']


    data = hdul_fwhm[1].data.field(0)
    hp = HEALPix(nside=nside, order=order, frame='icrs')

    pixsize = hp.pixel_resolution.to(u.arcsecond).value #steradian to square arcseconds  per pixel


    map_arc_s = ref_mapsize * ref_pixsize #map size in arcseconds
    npixside = ceil(map_arc_s / pixsize) #convert to map size in pixels for nu = 353 map.


    RA = np.linspace(center[0] - map_arc_s / 3600., center[0] + map_arc_s / 3600., npixside) * u.deg
    DEC = np.linspace(center[1] - map_arc_s / 3600. , center[1] + map_arc_s / 3600., npixside ) * u.deg

    #this creates the grid of RA / DEC values
    RA_grid, DEC_grid = np.meshgrid(RA, DEC)

    coords = SkyCoord(RA_grid.ravel(), DEC_grid.ravel(), frame='icrs')

    fwhm_arr = hp.interpolate_bilinear_skycoord(coords, data)

    fwhm = np.mean(fwhm_arr) #will want to use something better for later but this works for now..

    fwhm_pixels = fwhm * 60 / ref_pixsize #to get fwhm in pixels from arcmin since pixsize is arcseconds / pixel

    print(fwhm_pixels)
    fwhm_arcs = fwhm * 60

    print(fwhm_arcs)

    ###############Currently only using this to get the calfac conversion.
    bmsigma = fwhm_pixels / sqrt(8 * log(2))
    retext = round(fwhm_arcs * 5.0 / ref_pixsize)
    if retext % 2 == 0:
        retext += 1


    beam = Gaussian2DKernel(bmsigma , x_size=retext,y_size=retext, mode='oversample',factor=1)
    beam_norm = beam.array / np.sum(beam.array)

    calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0)))

    plt.imshow(beam_norm)
    plt.colorbar()
    plt.savefig('beam')
    plt.clf()

    to_MJy_Sr = 1 / (calfac * fwhm_arcs**2) #calibration factor based off of FWHM of our beam.
    return beam_norm, to_MJy_Sr

def read_in_fits(filename, center, ref_pixsize=8, ref_mapsize=260):

    hdul = fits.open(filename)
    data = hdul[1].data.field(0)
    head = hdul[1].header
    nside = head['NSIDE']
    order = head['ORDERING']



    hp = HEALPix(nside=nside, order=order, frame='icrs')
    #create a pixel grid in terms of the nu=353 grid for GNILC to create our intensity maps
    # sr = hp.pixel_area.to(u.degree)

    pixsize = hp.pixel_resolution.to(u.arcsecond).value


    map_arc_s = ref_mapsize * ref_pixsize #map size in arcseconds

    npixside = ceil(map_arc_s / pixsize) #convert to map size in pixels for nu = 353 map.

    RA = np.linspace(center[0] - map_arc_s / 3600., center[0] + map_arc_s / 3600., npixside) * u.deg
    DEC = np.linspace(center[1] - map_arc_s / 3600. , center[1] + map_arc_s / 3600., npixside ) * u.deg

    # this creates the grid of RA / DEC values
    RA_grid, DEC_grid = np.meshgrid(RA, DEC)

    coords = SkyCoord(RA_grid.ravel(), DEC_grid.ravel(), frame='icrs')

    # ipix = create_field_indexes([292.5, 89.95431463934247], ref_pixsize=8, ref_mapsize=260, nside=2048)
    map = hp.interpolate_bilinear_skycoord(coords, data)

    return map, pixsize

def blackbody_func(T, nu):
    h = 6.62 * 1e-34
    c = 3e8
    k = 1.38 * 1e-23
    wl = c / nu
    const =  h * c**2 / wl**5
    Tinv = np.divide(1, T)
    exponent = np.multiply(h * c / (k * wl), Tinv)
    denom = np.exp(exponent) - 1
    bb = np.divide(const, exponent)
    return bb

def calc_intensity(beta, tau,T, nu):
    nu_const = nu * 1e-9 / 353 #GHz
    power = np.power(nu_const, beta)
    bb = blackbody_func(T, nu)

    I = np.multiply(np.multiply(tau, bb), power)
    return I

def histograms(sim_data, real_data):

    sim = sim_data.flatten() - np.mean(sim_data)
    real = real_data.flatten() - np.mean(real_data)


    plt.hist(sim, int(50), label='sim', histtype='step')
    plt.hist(real, int(50), label='real', histtype='step')
    plt.xlabel('Flux [$\\frac{MJy}{Sr}$]')
    plt.ylabel('Nsources')
    plt.legend()
    plt.savefig('validation histograms')
    plt.clf()



if __name__ == '__main__':
    center = [206.9, -11.45]
    ref_size = 1000
    tau_name = 'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
    temp_name = 'COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
    beta_name = 'COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
    filenames = [tau_name, temp_name, beta_name]

    I_map =  create_map(filenames, ref_pixsize=8, ref_mapsize=ref_size, center=center)
    # convolved_I_map = convolve(I_map, beam, normalize_kernel=False)

    plt.title('Recreated Map')
    plt.imshow(I_map, origin='lower')
    plt.clim(0, 15)
    plt.colorbar().set_label('Flux [$\\frac{MJy}{sr}$]')
    plt.savefig('tester_boy.png')
    plt.clf()


    data, pixsize= read_in_fits('COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits', center, ref_pixsize=8, ref_mapsize=ref_size)
    side = int(sqrt(len(data)))
    map = np.reshape(data, (side, side))

    plt.title('Raw Data Map')
    plt.imshow(map, origin='lower')
    plt.clim(0, 15)
    plt.colorbar().set_label('Flux [$\\frac{MJy}{sr}$]')
    plt.savefig('tester_girl.png')
    plt.clf()

    histograms(I_map, map)
