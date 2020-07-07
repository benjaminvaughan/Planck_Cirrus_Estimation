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
from math_functions import *
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
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

def create_map(filenames, ref_head, ref_pixsize, ref_mapsize, center, nu):
    I_to_MJy = 1e20
    param_values = [] #Beta = 2, Temperature = 1, Tau = 0
    for f in filenames:
        data, pixsize, x_side, y_side, ra, dec = read_in_fits(f, center, ref_head, ref_pixsize, ref_mapsize)
        side = int(sqrt(len(data)))
        param_values.append(data)
    I = calc_intensity(param_values[2], param_values[0], param_values[1], nu)

    I_map = np.reshape(I, (x_side, y_side)) * I_to_MJy

    interped_map = interp_back_to_ref(I_map, ra, dec, ref_head, ref_mapsize)

    plt.imshow(interped_map)
    plt.savefig('../Test_Cases/test%scenter.png' % (center))
    plt.clf()

    return I_map
    # return interped_map

def read_in_fits(filename, center, ref_head, ref_pixsize=8, ref_mapsize=260):

    hdul = fits.open(filename)
    head = hdul[1].header
    if 'Temperature' in filename:
        data = hdul[1].data.field('TEMP')
        error = hdul[1].data.field('ERR_TEMP')
    elif 'Spectral-Index' in filename:
        data = hdul[1].data.field('BETA')
        error = hdul[1].data.field('ERR_BETA')

    elif 'Opacity' in filename:
        data = hdul[1].data.field('TAU353')
        error = hdul[1].data.field('ERR_TAU')

    else:
        data = hdul[1].data.field(0)
        print(head)
    nside = head['NSIDE']
    order = head['ORDERING']
    hdul.close()

    #Galactic Coordinate System
    hp = HEALPix(nside=nside, order=order, frame='galactic')
    #create a pixel grid in terms of the nu=353 grid for GNILC to create our intensity maps
    pixsize = hp.pixel_resolution.to(u.arcsecond).value

    map_arc_x = ref_mapsize[0] * ref_pixsize #map size in arcseconds
    map_arc_y = ref_mapsize[1] * ref_pixsize

    npixxside = ceil(map_arc_x / pixsize) #convert to map size in pixels for nu = 353 map.
    npixyside = ceil(map_arc_y / pixsize)


    x  = np.linspace(0, ref_mapsize[0],   npixxside)
    y  = np.linspace(0, ref_mapsize[1],   npixyside)

    X, Y = np.meshgrid(x, y)
    w = world(ref_head)
    skycoords = pixel_to_skycoord(X.ravel(), Y.ravel(), wcs=w, origin=0)
    RA_grid = np.asarray(skycoords.ra.to_string(decimal=True), dtype='float') * u.deg
    DEC_grid = np.asarray(skycoords.dec.to_string(decimal=True), dtype='float') * u.deg




    # coords = SkyCoord(RA_grid.ravel(), DEC_grid.ravel(), frame='icrs')
    coords = SkyCoord(ra=RA_grid.ravel(), dec=DEC_grid.ravel(), frame='icrs')
    gal_coords = coords.galactic

    map = hp.interpolate_bilinear_skycoord(gal_coords, data)

    x_side = len(x)
    y_side = len(y)
    return map, pixsize, y_side, x_side, RA_grid, DEC_grid


def interp_back_to_ref(img, ra, dec, ref_head, ref_shape):

    map_size = img.shape

    # reformat map data and coordinates
    data = np.ravel(img)
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_head)
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_shape[0], 0:ref_shape[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec))
    final_map = np.swapaxes(interp_map, 0, 1)
    return final_map


def get_fwhm(center, ref_pixsize=8, ref_mapsize=260):
    hdul_fwhm = fits.open('../Data/COM_CompMap_Dust-GNILC-Beam-FWHM_0128_R2.00.fits')

    header = hdul_fwhm[1].header
    nside = header['NSIDE']
    order = header['ORDERING']


    data = hdul_fwhm[1].data.field(0)
    hp = HEALPix(nside=nside, order=order, frame='galactic')

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
