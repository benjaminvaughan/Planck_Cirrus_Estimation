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

    return interped_map

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
    y = np.linspace(0, ref_mapsize[1],   npixyside)

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
    return interp_map
