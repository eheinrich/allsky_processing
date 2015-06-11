'''
equations_from_paper.py - Functions defining the projection equations used to
                          solve for the WCS.
                          
The projection equations used in this module were figured out by Borovicka et al., 
and can be found here: http://adsabs.harvard.edu/abs/1995A%26AS..112..173B

All coordinates must be expressed in radians.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
import numpy as np
import scipy

C = 1           # one of the fit parameters, but can be set=1 
                # to simplify the fitting process without losing information
                # its the global scale factor (shift of plate along the axis)

def r(x, y, params):
    '''
    Returns r as defined by Eq. 9
    
    '''
    a0 = params['a0']
    x0 = params['x0']
    y0 = params['y0']
    A  = params['A']
    F  = params['F']
    
    return C * (np.sqrt((x-x0)**2 + (y-y0)**2) + A * (y-y0) * np.cos(F-a0) - A * (x-x0) * np.sin(F-a0))
    
    
def u(r, params):
    '''
    Returns u as defined by Eq. 6    
    
    '''
    V  = params['V']
    S  = params['S']
    D  = params['D']
    P  = params['P']
    Q  = params['Q']

    return  V * r + S * (np.exp(D * r) - 1) + P * (np.exp(Q * r * r) - 1) 
    
    

def b(x, y, params):
    '''
    Returns b as defined by Eq. 4    
        
    '''
    a0 = params['a0']
    x0 = params['x0']
    y0 = params['y0']
    E  = params['E']
    
    return a0 - E + scipy.arctan2((y - y0),(x - x0))
    
    
def z(u, b, params):
    '''
    Return z as defined by Eq. 1
    
    '''
    e = params['e']        
    
    return np.arccos(np.cos(u) * np.cos(e) - np.sin(u) * np.sin(e) * np.cos(b))
    


def u_fit(star_x_y, a0, x0, y0, A, F, V, S, D, P, Q):
    '''
    Returns a u array for the given x,y array and parameter values.

    :params: numpy x,y array and projection parameters
    :return: Altitude numpy array

    '''
    params = {'a0' : a0,
              'x0' : x0,
              'y0' : y0,
              'A'  : A,
              'F'  : F,
              'V'  : V,
              'S'  : S,
              'D'  : D,
              'P'  : P,
              'Q'  : Q }

    r_fit = r(star_x_y[:,0], star_x_y[:,1], params)

    return u(r_fit, params)



def b_fit(star_x_y, a0, x0, y0):
    '''
    Returns an b array for the given x,y array and parameter values.
 
    :params: numpy x,y array and projection parameters
    :return: Azimuth numpy array     

    '''
    params = {'a0' : a0,
              'x0' : x0,
              'y0' : y0,
              'E'  : 0 }

    return b(star_x_y[:,0], star_x_y[:,1], params)
    

