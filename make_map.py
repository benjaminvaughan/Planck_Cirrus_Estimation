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

def makeGaussian(x_size, y_size, fwhm = 3, center=None):
    """ Make a gaussian kernel.

    size is the length of a side of the square.
    fwhm is full-width-half-maximum.
    center is where you want center of gaussian located, default is center of array.
    """

    sigma = fwhm / 2.355
    x = np.arange(0, x_size, 1, float)
    y = np.arange(0, y_size, 1, float)
    y = y[:,np.newaxis]

    if center is None:
        x0 = x_size // 2
        y0 = y_size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-1 * ((x-x0)**2 + (y-y0)**2) / (2*sigma**2))# / (2*np.pi*sigma**2)

def interpolator(img, ref_map, map):
    map_size = map['signal'].shape
    ref_size = ref_map['signal'].shape

    #create a grid of RA/DEC coordinates for image we want to interpolate
    w = world(map['shead'])

    #holder for the x, y pixel coordinates that we want.
    x = np.arange(0, map_size[0])
    y = np.arange(0, map_size[1])
    xvec = np.repeat(x[:, np.newaxis], map_size[1], axis=1)
    yvec = np.repeat(y[np.newaxis, :], map_size[0], axis=0)

    c = pixel_to_skycoord(xvec, yvec, w, origin=0)
    #this is converting the pixel coords to right ascension and declination in fk4
    ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
    dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)
    # reformat map data and coordinates
    data = np.ravel(img)
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_map['shead'])
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_size[0], 0:ref_size[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec), method='linear')
    return interp_map

def get_fwhm(filename, center, ref_pixsize=8, ref_mapsize=260):
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

    bmsigma = fwhm_pixels / sqrt(8 * log(2))
    retext = round(fwhm_arcs * 5.0 / ref_pixsize)
    if retext % 2 == 0:
        retext += 1

    beam = Gaussian2DKernel(bmsigma , x_size=retext,y_size=retext, mode='oversample',factor=1)
    beam_norm = beam.array / np.sum(beam.array)

    calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6)

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

def create_a_cut_out(RA, DEC):
    pass
    #need a way to calculate the cutout

def blackbody_func(T, nu):
    h = 6.62 * 1e-34
    c = 3e8
    k = 1.38 * 1e-23
    wl = c / nu
    const = 2 * h * c**2 / wl**5
    Tinv = np.divide(1, T)
    exponent = np.multiply(h * c / (k * wl), Tinv)
    denom = np.exp(exponent) - 1
    bb = np.divide(const, exponent)
    return bb

def calc_intensity(beta, tau,T, nu):
    # I_to_Jy = 1e-26
    nu_const = nu * 1e-9 / 353 #GHz
    power = np.power(nu_const, beta)
    bb = blackbody_func(T, nu)

    I = np.multiply(np.multiply(tau, bb), power) #* I_to_Jy
    return I

def histograms(sim_data, real_data):

    sim = sim_data.flatten() - np.mean(sim_data)
    real = real_data.flatten() - np.mean(real_data)


    plt.hist(sim, int(50), label='sim', histtype='step')
    plt.hist(real, int(50), label='real', histtype='step')
    plt.legend()
    plt.savefig('validation histograms')
    plt.clf()


if __name__ == '__main__':
    # center_ra = np.radians(360. - 206.9) #RIGHT ASCENSION IN DEGREES
    # center_dec = np.radians(-1 * -11.45 + 90) #Declination in degrees
    # center_ra = 10
    # center_dec = 120
    # center = [center_ra, center_dec]
    center = [206.9, -11.45]
    ref_size = 1000
    tau_name = 'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
    temp_name = 'COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
    beta_name = 'COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
    filename = [tau_name, temp_name, beta_name]
    param_values = [] #0 = Tau, 1 = Temp, 2 = Emissivity
    for f in filename:
        data, pixsize = read_in_fits(f, center, ref_pixsize=8, ref_mapsize=ref_size)
        side = int(sqrt(len(data)))
        param_values.append(data)
    nu = 353*1e9 #beta, tau,T, nu
    I = calc_intensity(param_values[2], param_values[0], param_values[1], nu)
    I_map = np.reshape(I, (side, side))

    beam, MJy_conv = get_fwhm(filename, center, ref_pixsize=pixsize, ref_mapsize=ref_size)

    convolved_I_map = convolve(I_map, beam, normalize_kernel=False)

    plt.title('Recreated Map')
    plt.imshow(convolved_I_map * MJy_conv * 1e6)
    plt.colorbar().set_label('$\\frac{MJy}{beam}$')
    plt.savefig('tester_boy.png')
    plt.clf()


    data, pixsize= read_in_fits('COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits', center, ref_pixsize=8, ref_mapsize=ref_size)
    side = int(sqrt(len(data)))
    map = np.reshape(data, (side, side))

    plt.title('Raw Data Map')
    plt.imshow(map)
    plt.colorbar().set_label('$\\frac{MJy}{sr}$ ')
    plt.savefig('tester_girl.png')
    plt.clf()

    histograms(convolved_I_map, map)
