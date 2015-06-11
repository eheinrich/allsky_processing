#!/usr/bin/env python
'''
calculation_steps.py - Functions to solve for the wcs projection parameters.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
from __future__ import division

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import logging
reload(logging)


import equations_from_paper as eq
from utils import astro_utils

logger = logging.getLogger(__name__)



def get_data(data_filename, ra_col, dec_col, x_col, y_col):
    '''
    Gets data of ra,dec and the corresponding x,y positions from file.
    
    :param data_filename: Filename containing this information
    :param ra_col: Column of data_filename containing ra
    :param dec_col: Column of data_filename containing dec
    :param x_col: Column of data_filename containing x positions
    :param y_col: Column of data_filename containing y positions
    :return star_ra_dec: numpy array with ra, dec vals
    :return star_x_y: numpy array with x, y vals

    '''
    star_ra_dec = np.genfromtxt(data_filename, usecols = (ra_col, dec_col))
    star_x_y    = np.genfromtxt(data_filename, usecols = (x_col, y_col))

    return star_ra_dec, star_x_y
    


def az_zendist_of_data(image_filename, ra_dec, site=''):
    '''
    Returns an array with azimuth/zenith distance locations corresponding to
    ra,dec information.

    :param image_filename: Path of FITS file of image
    :param ra_dec: Array of ra, dec information [decimal degrees]
    :return: Array of azimuth, zenith distance in [radians]

    '''
    obstime = fits.getval(image_filename, 'DATE-OBS').replace('T', ' ')
    
    if not site:
        site = fits.getval(image_filename, 'SITE')

    observer = astro_utils.observer(site, obstime)
    
    az_zendist = np.zeros(shape=ra_dec.shape)
    
    for i in range(len(az_zendist)):
        az_zen_am = np.deg2rad(astro_utils.world2az_zendist(ra_dec[i,0], ra_dec[i,1], observer))
        az_zendist[i,0] = az_zen_am[0]
        az_zendist[i,1] = az_zen_am[1]

    return az_zendist


  
def initial_fits(x_y, az_zendist, p0):
    '''
    Does initial fits to get better guesses on the parameters.
    
    :param x_y: Array containing x,y pixel positions
    :param az_zendist: Array with azimuth, zenith distance information [radians]
    :param p0: Input dictionary of parameters
    :return: Output dictionary of parameters
    
    '''    
    logger.info('Doing initial fits to determine better guesses for a0, x0, y0, A, F, and the lens constants.')
    logger.info('Starting with params: ')
    logger.info(p0)
    
    logger.info('Commencing initial azimuth fits')
    p0_az = (p0['a0'], p0['x0'], p0['y0'])    
    fit_az, cov_az = curve_fit(eq.b_fit, x_y, az_zendist[:,0], p0=p0_az) 
    logger.debug('first initial fit: ', fit_az)
    
    logger.info('Commencing initial zenith distance fits')
    p0_alt = (fit_az[0], fit_az[1], fit_az[2], p0['A'], p0['F'], p0['V'], p0['S'], p0['D'], p0['P'], p0['Q'])
    fit_alt, cov_alt = curve_fit(eq.u_fit, x_y, az_zendist[:,1], p0=p0_alt)
    logger.debug('second_initial fit: ', fit_alt)
    
    # reassign initial parameters
    p0['a0'] = fit_alt[0]
    p0['x0'] = fit_alt[1]
    p0['y0'] = fit_alt[2]
    p0['A']  = fit_alt[3]
    p0['F']  = fit_alt[4]
    p0['V']  = fit_alt[5]
    p0['S']  = fit_alt[6]
    p0['D']  = fit_alt[7]
    p0['P']  = fit_alt[8]
    p0['Q']  = fit_alt[9]
    
    logger.info('Finished initial fits with params: ')
    logger.info(p0)
    
    return p0    



def guess_for_E_x0_y0(x_y, az_zendist, params):
    '''
    Finds a guess for E. After E is refined, x0 and y0 are recalculated.

    :param x_y: Array of x,y pixel positions
    :param az_zendist: Array of azimuth, zenith distance [radians]
    :param p0: Input dictionary of parameters
    :return: Output dictionary of parameters
    
    '''
    logger.info('Improving guesses for E, x0, y0')
    
    E_tries     = np.linspace(-4*np.pi, 4*np.pi, 500)
    E_tries_min = np.zeros(shape=E_tries.shape)     
    numTries = E_tries.shape[0]

    x0_old, y0_old = params['x0'], params['y0']

    for i in range(numTries):
        params['x0'] = x0_old + (params['e'] / params['V']) * np.cos(E_tries[i])
        params['y0'] = y0_old + (params['e'] / params['V']) * np.sin(E_tries[i])
        params['E']  = E_tries[i]

        E_tries_min[i] = minimize_initial_fits(x_y, az_zendist, params) 
            
    E_min = E_tries[np.argmin(E_tries_min)]

    params['E']  = E_min
    params['x0'] = x0_old + (params['e'] / params['V']) * np.cos(E_min)
    params['y0'] = y0_old + (params['e'] / params['V']) * np.sin(E_min)

    logger.info('Final best guess of all parameters: ')
    logger.info(params)

    return params



def calculateAz(arg, answ, a, E):
    '''
    Calculates the azimuth.    
    
    :param arg: Argument to arcsin of Eq. 2 -- sin(b)*sin(u)/sin(z)
    :param answ: Absolute val of the arcsin of arg -- abs(arcsin(arg))
    :param a: Measured azimuth
    :param E: Projection parameter E
    :return: Calculated azimuth

    '''  
    if arg < 0:
        a1 = np.pi + answ + E
        a2 = 2*np.pi - answ + E     
    elif arg > 0:
        a1 = answ + E
        a2 = np.pi - answ + E             
    elif arg == 0:
        a1 = E
        a2 = np.pi + E
    else:
        raise ValueError('Unexpected value for arg, got %s. Maybe try decreasing the step size (fraction).' % str(arg))
        
    da1 = np.abs(a - a1)
    da2 = np.abs(a - a2)
    
    a_diff_min = min(da1, da2)
    
    if a_diff_min == da1:
        a_calc = a1
    if a_diff_min == da2:
        a_calc = a2
            
    return a_calc
    
    
    
def find_azimuth_zendist(x, y, a_actual, params):
    '''
    Returns calculated azimuth and altitude, as well as u, b (altitude and azimuth
    of the camera axis)  
    
    :param x: x pixel position
    :param y: y pixel position
    :param a_actual: Measured azimuth
    :param p: Dictionary of projection parameters
    :return: Calculated azimuth, altitude, b, and u

    '''
    r = eq.r(x, y, params)
    u = eq.u(r, params)
    b = eq.b(x, y, params)
    z = eq.z(u, b, params)
    
    # solve for azimuth
    arg  = np.sin(b) * np.sin(u) / np.sin(z)
    answ = abs(np.arcsin(arg))
    a    = calculateAz(arg, answ, a_actual, params['E'])

    return (a, z, b, u)



def calc_az_zendist_with_params(x_y, az_zendist, params):
    '''
    Calculate azimuth, zenith distance given a set of parameters.
    
    :param x_y: x,y pixel numpy array
    :param az_zendist: real azimuth, zenith distance numpy array
    :param params: dictionary of projection parameters
    :return: calculated azimuth, zenith distance postions based on the given parameters,
             as well as the residuals of the positions
    
    '''
    calc_az_zendist = np.zeros(shape=az_zendist.shape)    
    for i in range(len(az_zendist)):
        calc_az_zendist[i,0] , calc_az_zendist[i,1], w, w = find_azimuth_zendist(x_y[i,0], x_y[i,1], az_zendist[i,0], params) 
        
    residuals = (az_zendist - calc_az_zendist)
    
    return calc_az_zendist, residuals



def minimize_initial_fits(x_y, az_zendist, params):
    '''
    Function to minimize for initial fits.    
    
    :params: x, y, azimuth, and projection parameters
    :return: Sum to be minimized

    '''
    calc_az_zendist, residuals = calc_az_zendist_with_params(x_y, az_zendist, params)

    F_G = np.subtract(az_zendist, calc_az_zendist)    

    F_G_abs  = np.abs(F_G)
    finalSum = np.sum(F_G_abs)
    
    return finalSum
    
    

def interactive_3D_plot(x_y, az_zendist, calc_az_zendist):
    '''
    Makes a 3D matplotlib plot, comparing measured alt,az and calculated alt,az
    
    :param x_y: Array of x,y positions
    :param azimuth_zendist: Array of measured azimuth,zenith distance
    :param calc_az_alt: array of calculated az,alt
    
    '''
    #Plot alt
    fig1 = plt.figure()

    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.scatter(x_y[:,0], x_y[:,1], az_zendist[:,1], c='r')    
    ax1.scatter(x_y[:,0], x_y[:,1], calc_az_zendist[:,1], c = 'b')    

    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Altitude')
    
    #Plot az
    fig2 = plt.figure()
    
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(x_y[:,0], x_y[:,1], az_zendist[:,0], c='r')
    ax2.scatter(x_y[:,0], x_y[:,1], calc_az_zendist[:,0], c = 'b')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Azimuth')
    
    plt.show()
    
    return
    
    
    