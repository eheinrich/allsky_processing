'''
solved_parameters.py - Contains the solved wcs projection parameters for each site.

The projection parameters aren't expected to change unless the camera is moved. Also note that 
they may appear to be incorrect if the time stamps on the images is incorrect.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
from math import radians, pi

mm2pix = 0.0074 # conversion factor

_sites = {'LSC' : ['-70.815',        '-30.165',       2215],
          'NGQ' : ['80.0166667',     '32.31666',      5100],
          'OGG' : ['-156.258055556', '20.7069444444', 3065]}



def parameters(site):
    '''
    Returns WCS projection parameters for a given site.

    :param site: str, can be 'LSC', 'OGG', 'NGQ', etc. 
                 Setting site = 'generic', will return initial guesses of the parameters
                 to be used to solve for them.
                 
    :return: [a0_i, x0_i, y0_i, A_i, F_i, V_i, S_i, D_i, P_i, Q_i, e, E], 
             [FR, Rtod, theta(in radians), x0_i, y0_i], 
             [str longitude(- for west, + for east, degrees), str latitude(- for south, + for north, degfrees), float elevation(meters)]

    '''  
    if site == 'generic':                           # resonable initial guesses on projection parameters
        parameters = {'a0'  : -1,                   # angle between the x axis and the direction to the south
                      'x0'  : 315,                  # center of x-coord projection on xy plane
                      'y0'  : 250,                  # center of y-coord projection on xy plane
                      'A'   : 0,                    # related to eccentricity of proj due to tilt of camera
                      'F'   : -15,                  # related to direction of inclination of projection due to tilt of camera
                      'V'   : 0.72 * mm2pix,        # lens constant
                      'S'   : 0.0009,               # lens constant
                      'D'   : 0.6 * mm2pix,         # lens constant
                      'P'   : 0,                    # lens constant
                      'Q'   : 0,                    # lens constant
                      'e'     : 0.002,
                      'E'     : 0.1,
                      'theta' : radians(-2.4),    
                      'f_rad' : 302.5/90,
                      'r2d'   : 90/(pi/2)}   
                      
    if site == 'OGG':                                 # hawaii
        parameters = {'a0'    : -1.52839897e+00,
                      'x0'    : 3.29605137e+02,
                      'y0'    : 2.32098951e+02,
                      'A'     : 8.54863794e-04,
                      'F'     : -1.63511077e+01,
                      'V'     : 5.20268597e-03,
                      'S'     : -2.95160214e-04,
                      'D'     : -3.30385847e+01,
                      'P'     : -7.95110819e+05,
                      'Q'     : 1.17331349e-13,
                      'e'     : 1.09127438e-03,
                      'E'     : 9.43083039e-04,
                      'theta' : radians(-2.4),    
                      'f_rad' : 302.5/90,
                      'r2d'   : 90/(pi/2),
                      'lon'   : '-156.258055556',
                      'lat'   : '20.7069444444',
                      'elev'  : 3065.0}
	
    if site == 'LSC':                                 # chile
        parameters = {'a0_i'  : -1.28125108e+00,
                      'x0_i'  : 3.14744253e+02,
                      'y0_i'  : 2.49567438e+02,
                      'A_i'   : -2.54309250e-03,
                      'F_i'   : -1.59032966e+01,
                      'V_i'   :  5.16687135e-03,
                      'S_i'   : -9.41172090e-04,
                      'D_i'   : -3.30385847e+01,
                      'P_i'   : -7.95110819e+05,
                      'Q_i'   : 4.83936764e-14,
                      'e'     : -1.52825089e-02,
                      'E'     : 6.23623074e-02,
                      'theta' : radians(-9.7),    
                      'f_rad' : 302.5/90,
                      'r2d'   : 90/(pi/2),
                      'lon'   : '-70.804672222',
                      'lat'   : '-30.1673472222',
                      'elev'  : 2198.0}
                      
    if site == 'NGQ':                                 # china
        parameters = {'a0'    : -3.90086629,          # projection parameters updated on 2015-06-11 EHJ
                      'x0'    : 316.13044633,         
                      'y0'    : 242.2933928,
                      'A'     : 0.00114328,
                      'F'     : 12.70014376,
                      'V'     : 0.00519457,
                      'S'     : 0.75639737,
                      'D'     : -5.62003607,
                      'P'     : -90572.87223312,
                      'Q'     : -64.28254367,
                      'e'     : -0.01822852,
                      'E'     : -0.00097206,
                      'theta' : radians(-226.7),    
                      'f_rad' : 302.5/90,
                      'r2d'   : 90/(pi/2),
                      'lon'   : '80.0166667',
                      'lat'   : '32.31666',
                      'elev'  : 5100.0}

    if not parameters:
        raise NameError('WCS fit parameters for site %s are not defined' % site)
        return 1

    return parameters


