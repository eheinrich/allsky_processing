'''
astro_utils.py - Functions to perform various astro calculations.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
from __future__ import division

import ephem
import numpy as np

from wcs import solved_parameters



def observer(site, obstime):
    '''
    Returns a pyephem observer object.
    
    :param site:    Three letter site code, lowercase
    :param obstime: Date of the observer, i.e. '2015-05-29 20:05:58.359'
    
    :return: Observer object
    '''
    observer = ephem.Observer()
    observer.lon, observer.lat, observer.elevation = solved_parameters._sites[site]
    observer.date = obstime

    return observer



def world2az_zendist(ra, dec, observer):
    '''
    Returns [azimuth, altitude, airmass] for a given ra, dec, observer.
    
    :param ra:       Right ascension J2000, degrees
    :param dec:      Declinations J2000, degrees
    :param observer: Pyephem observer object
    
    :return:         Corresponding azimuth [deg], altitude [deg], airmass
    
    '''
    star = ephem.FixedBody()
    
    star._ra  = str(ra)
    star._dec = str(dec)
    star.compute(observer)

    return [np.degrees(star.az), 
            90 - np.degrees(star.alt), 
            airmass(90-np.degrees(star.alt))]



def flux_to_mag(flux, dflux):
    '''
    Returns the instrumental magnitude and its error given flux and flux error.
    The flux should be in counts per second.

    '''   
    mag  = -2.5 * np.log10(flux)
    dmag =  2.5 * dflux / (flux * np.log(10))

    return [mag, dmag]



def airmass(zen_dist):
    '''
    Returns airmass for apparent zenith distance with truncation at > 85 deg
    zenith. This particular airmass formula can be found on the Wikipedia
    page on airmasses.
    
    :param zen_dist:  Zenith distance, degrees
    
    :return:          Airmass
    
    '''
    zen_dist = abs(np.radians(zen_dist)) if zen_dist <= 85 else np.radians(85)
    
    secant_z = 1/np.cos(zen_dist)
        
    return secant_z - (0.0018167 * (secant_z-1) + 0.002875 * (secant_z-1)**2 + 0.0008083 * (secant_z-1)**3)



def moonTracker(observer):
    '''
    Returns az, alt of the Moon for a given observer. 
    
    :param observer: Pyephem observer object
    :return: Moon azimuth, altitude [degrees]
    '''
    moon = ephem.Moon()
    moon.compute(observer)    
    
    return (np.degrees(moon.az), np.degrees(moon.alt))

    
    
def visiblePlanets(observer):
    '''
    Returns major planets above the horizon for a give observer.
    In order: Jupiter, Mars, Saturn, Venus, Mercury.
    List returned gives for each planet the name, azimuth, altitude.

    :param observer: Pyephem observer object
    :return: list of planet positions [deg] and names, i.e. list of ['Jupiter', azimuth, altitude] ...
    
    '''
    jupiter = ephem.Jupiter()
    mars = ephem.Mars()
    saturn = ephem.Saturn()
    venus = ephem.Venus()
    mercury = ephem.Mercury()
    
    jupiter.compute(observer)
    mars.compute(observer)
    saturn.compute(observer)
    venus.compute(observer)
    mercury.compute(observer)
    
    planetList = [['Jupiter', np.degrees(jupiter.az), np.degrees(jupiter.alt)], 
                  ['Mars',    np.degrees(mars.az),    np.degrees(mars.alt)   ], 
                  ['Saturn',  np.degrees(saturn.az),  np.degrees(saturn.alt) ],
                  ['Venus',   np.degrees(venus.az),   np.degrees(venus.alt)  ],
                  ['Mercury', np.degrees(mercury.az), np.degrees(mercury.alt)]] 
    return planetList



def dist_to_pix(distance, pixscale, u='degree'):
    '''
    Converts to distance in pixels    
    
    :param distance: distance in radians
    :param pixscale: degrees per pixel
    :param u: unit of distance, either 'radian' or 'degree'
    :return: distance in pixels
    
    '''
    if u == 'radian':
        pixdist = np.degrees(distance) / pixscale
    if u == 'degree':
        pixdist = distance / pixscale

    return pixdist
    
    
    