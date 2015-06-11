'''
derivatives.py - Contains the derivatives and calculations necessary to do the least-squares 
                 fitting to solve for the WCS projection.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
from __future__ import division

import logging
reload(logging)

from numpy import sin, cos, negative, exp, sqrt, zeros, ones, linalg

from calc_steps import find_azimuth_zendist
import equations_from_paper as eq


logger = logging.getLogger(__name__)



class Derivatives(object):
    '''
    Contains all the derivatives and the loop that linearizes and solves for the 
    parameters.
    
    '''
    C = 1
    
    def __init__(self, params):
        # initialize the parameters
        self.p        = params
        
        # initialize the following with dummy values
        # will make the methods look cleaner
        self.xpix   = 1
        self.ypix   = 1
        self.r      = 1
        self.a_real = 1
        self.a_calc = 1
        self.z_real = 1
        self.z_calc = 1
        self.b      = 1
        self.u      = 1


    def update_param_vals(self, dp, fraction):
        self.p['a0'] += fraction*dp[0]
        self.p['x0'] += fraction*dp[1]
        self.p['y0'] += fraction*dp[2]
        self.p['A']  += fraction*dp[3]
        self.p['F']  += fraction*dp[4]
        self.p['V']  += fraction*dp[5]
        self.p['S']  += fraction*dp[6]
        self.p['D']  += fraction*dp[7]
        self.p['P']  += fraction*dp[8]
        self.p['Q']  += fraction*dp[9]
        self.p['e']  += fraction*dp[10]
        self.p['E']  += fraction*dp[11]


    def leastSqSolve(self, x_y, az_zendist, fraction, param_tolerance):
        '''
        Does the least squares fitting. 
        
        :param x_y: x,y pixel positions
        :param az_zendist: Azimuth, zenith distance positions
        :param fraction: float multiplier to determine step size
        :return: Fitted parameters, uncertainty of parameters
        
        '''
        logger.info('Beginning least squares solve of parameters.')    
    
        numSources, numParams = len(az_zendist), len(self.p)
        
        # initialize arrays used in the fitting process
        tolerance_dp = negative(ones(shape=(numParams, 1)))
        f_g_array    = zeros(shape=(2*numSources, 1))
        dp           = zeros(shape=(numParams, 1)) 
        derivMatrix  = zeros(shape=(2*numSources, numParams))
        
        while (abs(dp) > tolerance_dp).any():
            try: 
                tolerance_dp = param_tolerance * ones(shape=(numParams, 1))     
    
                logger.info('New loop parameters, parameter errors: ')
                logger.info('a0   %10.7f\t%3.7f' % (self.p['a0'], dp[0]))
                logger.info('x0   %10.7f\t%3.7f' % (self.p['x0'], dp[1]))
                logger.info('y0   %10.7f\t%3.7f' % (self.p['y0'], dp[2]))
                logger.info('A    %10.7f\t%3.7f' % (self.p['A'],  dp[3]))
                logger.info('F    %10.7f\t%3.7f' % (self.p['F'],  dp[4]))
                logger.info('V    %10.7f\t%3.7f' % (self.p['V'],  dp[5]))
                logger.info('S    %10.7f\t%3.7f' % (self.p['S'],  dp[6]))
                logger.info('D    %10.7f\t%3.7f' % (self.p['D'],  dp[7]))
                logger.info('P    %10.7f\t%3.7f' % (self.p['P'],  dp[8]))
                logger.info('Q    %10.7f\t%3.7f' % (self.p['Q'],  dp[9]))
                logger.info('e    %10.7f\t%3.7f' % (self.p['e'],  dp[10]))
                logger.info('E    %10.7f\t%3.7f\n' % (self.p['E'],  dp[11]))
    
                self.update_param_vals(dp, fraction)
                                                
                for i in range(numSources):
                
                    self.xpix   = x_y[i,0]
                    self.ypix   = x_y[i,1]            
                    self.a_real = az_zendist[i,0]
                    self.z_real = az_zendist[i,1]
        
                    self.a_calc, self.z_calc, self.b, self.u = find_azimuth_zendist(self.xpix, self.ypix, self.a_real, self.p)
                    self.r = eq.r(self.xpix, self.ypix, self.p)
    
                    derivMatrix[i,0]  = self.df_da0()
                    derivMatrix[i,1]  = self.df_dx0()
                    derivMatrix[i,2]  = self.df_dy0()
                    derivMatrix[i,3]  = self.df_dA()
                    derivMatrix[i,4]  = self.df_dF()
                    derivMatrix[i,5]  = self.df_dV()
                    derivMatrix[i,6]  = self.df_dS()
                    derivMatrix[i,7]  = self.df_dD()
                    derivMatrix[i,8]  = self.df_dP()
                    derivMatrix[i,9]  = self.df_dQ()
                    derivMatrix[i,10] = self.df_de()
                    derivMatrix[i,11] = self.df_dE()
                    
                    f_g_array[i] = -self.f()  
                    
                for i in range(numSources, 2*numSources):
                
                    self.xpix   = x_y[i-numSources,0]
                    self.ypix   = x_y[i-numSources,1]
                    self.a_real = az_zendist[i-numSources,0]
                    self.z_real = az_zendist[i-numSources,1]
                    
                    self.a_calc, self.z_calc, self.b, self.u = find_azimuth_zendist(self.xpix, self.ypix, self.a_real, self.p)
                    self.r = eq.r(self.xpix, self.ypix, self.p)
    
                    derivMatrix[i,0]  = self.dg_da0()
                    derivMatrix[i,1]  = self.dg_dx0()
                    derivMatrix[i,2]  = self.dg_dy0()
                    derivMatrix[i,3]  = self.dg_dA()
                    derivMatrix[i,4]  = self.dg_dF()
                    derivMatrix[i,5]  = self.dg_dV()
                    derivMatrix[i,6]  = self.dg_dS()
                    derivMatrix[i,7]  = self.dg_dD()
                    derivMatrix[i,8]  = self.dg_dP()
                    derivMatrix[i,9]  = self.dg_dQ()
                    derivMatrix[i,10] = self.dg_de()
                    derivMatrix[i,11] = self.dg_dE()
                    
                    f_g_array[i] = -self.g()
                
                dp, residuals, rank, singular_vals = linalg.lstsq(derivMatrix, f_g_array)

            except KeyboardInterrupt:
                return self.p, dp

        return self.p, dp
        

    def f(self):
        return abs(self.a_calc - self.a_real) * sin(self.z_real)
                
    def g(self):
        return abs(self.z_calc - self.z_real)
                
    def df_da(self):
        if self.a_calc > self.a_real:
            deriv = sin(self.z_real)
        if self.a_calc < self.a_real:
            deriv = -sin(self.z_real)
        else:
            deriv = sin(self.z_real)
    
        return deriv
                
    def dg_dz(self):
        if self.z_calc > self.z_real:
            deriv = 1.
        if self.z_calc < self.z_real:
            deriv = -1.
        else:
            deriv = 1.
    
        return deriv
            
    def da_db(self):
        return (cos(self.b) * sin(self.u)) / (cos(self.a_calc - self.p['E']) * sin(self.z_calc))
        
    def da_dz(self):
        return negative(sin(self.b) * sin(self.u) * cos(self.z_calc)) / (cos(self.a_calc - self.p['E']) * sin(self.z_calc) * sin(self.z_calc))
        
    def da_du(self):
        return (sin(self.b) * cos(self.u)) / (cos(self.a_calc - self.p['E']) * sin(self.z_calc))
        
    def da_dE(self):
        return 1
        
    def dz_db(self):
        return negative(sin(self.u) * sin(self.p['e']) * sin(self.b)) / sin(self.z_calc)    
        
    def dz_du(self):
        return (sin(self.u) * cos(self.p['e']) + cos(self.u) * sin(self.p['e']) * cos(self.b)) / sin(self.z_calc)    
        
    def dz_de(self):
        return (cos(self.u) * sin(self.p['e']) + sin(self.u) * cos(self.p['e']) * cos(self.b)) / sin(self.z_calc)    
        
    def db_da0(self):
        return 1    
        
    def db_dx0(self):
        y_diff = self.ypix - self.p['y0']
        x_diff_squared = (self.xpix - self.p['x0'])**2
        return (y_diff / x_diff_squared) * (cos(self.b - self.p['a0'] + self.p['E']))**2
            
    def db_dy0(self):
        x_diff = self.xpix - self.p['x0']
        return negative((cos(self.b - self.p['a0'] + self.p['E']))**2) / x_diff
            
    def db_dE(self):
        return -1
            
    def du_dr(self):
        l = self.p['V']
        m = self.p['S'] * self.p['D'] * exp(self.p['D'] * self.r)
        n = 2 * self.p['Q'] * self.p['P'] * self.r * exp(self.p['Q'] * self.r * self.r)
        return l + m + n    
        
    def du_dV(self):
        return self.r
            
    def du_dS(self):
        return exp(self.p['D'] * self.r) - 1
    
    def du_dD(self):
        return self.p['S'] * self.r * exp(self.p['D'] * self.r)
            
    def du_dP(self):
        return exp(self.p['Q'] * self.r * self.r) - 1    
        
    def du_dQ(self):
        return self.p['P'] * self.r * self.r * exp(self.p['Q'] * self.r * self.r)    
        
    def dr_dx0(self):
        x_diff = self.xpix - self.p['x0']
        y_diff = self.ypix - self.p['y0']
        r_simple = sqrt(x_diff**2 + y_diff**2)    
        return self.C * ((-x_diff / r_simple) + self.p['A'] * sin(self.p['F'] - self.p['a0']))
    
    def dr_dy0(self):
        x_diff = self.xpix - self.p['x0']
        y_diff = self.ypix - self.p['y0']
        r_simple = sqrt(x_diff**2 + y_diff**2) 
        return self.C * ((-y_diff / r_simple) - self.p['A'] * cos(self.p['F'] - self.p['a0']))
        
    def dr_da0(self):
        x_diff = self.xpix - self.p['x0']
        y_diff = self.ypix - self.p['y0']
        return self.p['A'] * self.C * (y_diff * sin(self.p['F'] - self.p['a0']) + x_diff * cos(self.p['F'] - self.p['a0']))
        
    def dr_dF(self):    
        x_diff = self.xpix - self.p['x0']
        y_diff = self.ypix - self.p['y0']
        return -self.p['A'] * self.C * (y_diff * sin(self.p['F'] - self.p['a0']) + x_diff * cos(self.p['F'] - self.p['a0']))
        
    def dr_dA(self):
        x_diff = self.xpix - self.p['x0']
        y_diff = self.ypix - self.p['y0']
        return self.C * (y_diff * cos(self.p['F'] - self.p['a0']) - x_diff * sin(self.p['F'] - self.p['a0']))
                
    def df_da0(self):
        l = self.da_db() * self.db_da0()
        m = self.da_du() * self.du_dr() * self.dr_da0()
        n = self.da_dz() * (self.dz_du() * self.du_dr() * self.dr_da0() + self.dz_db() * self.db_da0())
        return self.df_da() * (l + m + n)    
        
    def dg_da0(self):
        a = self.dz_du() * self.du_dr() * self.dr_da0()
        b = self.dz_db() * self.db_da0()
        return self.dg_dz() * (a + b)
            
    def df_dx0(self):
        l = self.da_db() * self.db_dx0()
        m = self.da_du() * self.du_dr() * self.dr_dx0()
        n = self.da_dz() * (self.dz_du() * self.du_dr() * self.dr_dx0() + self.dz_db() * self.db_dx0())
        return self.df_da() * (l + m + n)    
        
    def df_dy0(self):
        l = self.da_db() * self.db_dy0()
        m = self.da_du() * self.du_dr() * self.dr_dy0()
        n = self.da_dz() * (self.dz_du() * self.du_dr() * self.dr_dy0() + self.dz_db() * self.db_dy0())
        return self.df_da() * (l + m + n)
            
    def dg_dx0(self):
        l = self.dz_du() * self.du_dr() * self.dr_dx0()
        m = self.dz_db() * self.db_dx0()
        return self.dg_dz() * (l + m)
            
    def dg_dy0(self):
        l = self.dz_du() * self.du_dr() * self.dr_dy0()
        m = self.dz_db() * self.db_dy0()
        return self.dg_dz() * (l + m)    
    
    def df_dV(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m)
        return base * self.du_dV()
            
    def dg_dV(self):
        base = self.dg_dz() * self.dz_du()
        return base * self.du_dV()
            
    def df_dS(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m)
        return base * self.du_dS()
            
    def dg_dS(self):
        base = self.dg_dz() * self.dz_du()
        return base * self.du_dS()
        
    def df_dD(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m)
        return base * self.du_dD()
            
    def dg_dD(self):
        base = self.dg_dz() * self.dz_du() * self.du_dr()
        return base * self.du_dD()
        
    def df_dP(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m)
        return base * self.du_dP()         
    
    def dg_dP(self):
        base = self.dg_dz() * self.dz_du()
        return base * self.du_dP()    
    
    def df_dQ(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m)
        return base * self.du_dQ()
    
    def dg_dQ(self):
        base = self.dg_dz() * self.dz_du()
        return base * self.du_dQ()
            
    def df_de(self):
        return self.df_da() * self.da_dz() * self.dz_de()
            
    def dg_de(self):
        return self.dg_dz() * self.dz_de()
            
    def df_dE(self):
        l = self.da_dE()
        m = (self.da_db() + self.da_dz() * self.dz_db()) * self.db_dE()
        return self.df_da() * (l + m)
            
    def dg_dE(self):
        return self.dg_dz() * self.dz_db() * self.db_dE()
            
    def df_dA(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m) * self.du_dr()
        return base * self.dr_dA()
            
    def dg_dA(self):
        base = self.dg_dz() * self.dz_du() * self.du_dr()
        return base * self.dr_dA()
            
    def df_dF(self):
        l = self.da_du()
        m = self.da_dz() * self.dz_du()
        base = self.df_da() * (l + m) * self.du_dr()
        return base * self.dr_dF()
            
    def dg_dF(self):
        base = self.dg_dz() * self.dz_du() * self.du_dr()
        return base * self.dr_dF()
    
    
    