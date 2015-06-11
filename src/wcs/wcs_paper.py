#!/opt/anaconda/bin/python

import numpy as np
import ephem as ep
import scipy.optimize as op
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy as sc
from astropy.io import fits

def main():
    image_data = '/home/eheinrich/Desktop/Code/CloudsCode/china/Inputs/stars_china.dat'
    current_image = '/home/eheinrich/Desktop/Code/CloudsCode/china/Inputs/allsky_china.fits'
    
    #SITE INFORMATION
#    china_obs_lon, china_obs_lat, china_obs_elev = '97.73055556', '37.37777778', 3200
    china_obs_lon, china_obs_lat, china_obs_elev = '80.0166667', '32.31666', 5100

#    chile_obs_lon, chile_obs_lat, chile_obs_elev = '-70.815', '-30.165', 2215
    
    #GET OBSERVATION TIME FROM HEADER OF IMAGE USED
    hdulist = fits.open(current_image)
    OBSTIME = hdulist[0].header['DATE-OBS']
    hdulist.close()
    OBSTIME = OBSTIME.replace('T', ' ')
    
    #specify location of observation
    obs = ep.Observer()
    obs.lon, obs.lat, obs.elevation = china_obs_lon, china_obs_lat, china_obs_elev
    obs.date = OBSTIME
    
    mm2pix = 0.0074
    C = 1.
    
    ##FROM A STARS.DAT-LIKE FILE
    star_ra_dec = np.genfromtxt(image_data, usecols = (1,2))
    star_x_y = np.genfromtxt(image_data, usecols = (3,4))
    ##star_names = np.genfromtxt('stars.dat', usecols = 0, dtype = str )
    star_az_alt = np.zeros(shape = star_ra_dec.shape)
    
    #FROM A MATCHED_STARS.TXT-LIKE FILE
    #star_ra_dec = np.genfromtxt(image_data, usecols=(1,2))
    #star_x_y = np.genfromtxt(image_data, usecols=(6,8))
    #star_names = np.genfromtxt('stars.dat', usecols = 0, dtype = str )
    #star_az_alt = np.zeros(shape = star_ra_dec.shape)
    
    #number of sources
    numSources = star_ra_dec.shape[0]
    
    #compute az/alt
    for i in range(numSources):    
        star = ep.FixedBody()
        star._ra = np.str(star_ra_dec[i,0])
        star._dec = np.str(star_ra_dec[i,1])
        star.compute(obs)
        star_az_alt[i,0] = np.deg2rad(np.degrees(star.az))
        star_az_alt[i,1] = np.pi/2 - np.deg2rad(np.degrees(star.alt))
    
    #Define functions   
    def alt_initial_fit(star_x_y, a0, x0, y0, A, F, V, S, D, P, Q):
        r = C * (np.sqrt(((star_x_y[:,0]-x0)**2) + ((star_x_y[:,1]-y0)**2)) + A * (star_x_y[:,1]-y0) * np.cos(F - a0) - A * (star_x_y[:,0] - x0) * np.sin(F - a0))
        return V * r + S * (np.exp(D * r) - 1.) + P * (np.exp(Q * r * r) - 1.)
        
    def az_initial_fit(star_x_y, a0, x0, y0):
        return a0 + sc.arctan2((star_x_y[:,1] - y0),(star_x_y[:,0] - x0)) #+ np.pi
    
    def find_az(E, a):
        return np.sin(a - E)  
        
    def determineAz(arg, answ, a_actual, E):
        
        if arg < 0.:
            a1 = np.pi + answ + E
            a2 = 2 * np.pi - answ + E
            
            a1diff = np.abs(a_actual - a1)
            a2diff = np.abs(a_actual - a2)
            a_diff_min = min(a1diff, a2diff)
            
            if a_diff_min == a1diff:
                a_calc = a1
            if a_diff_min == a2diff:
                a_calc = a2
                
        if arg > 0.:
            a1 = answ + E
            a2 = np.pi - answ + E
            
            a1diff = np.abs(a_actual - a1)
            a2diff = np.abs(a_actual - a2)
            a_diff_min = min(a1diff, a2diff)
            
            if a_diff_min == a1diff:
                a_calc = a1
            if a_diff_min == a2diff:
                a_calc = a2
                
        if arg == 0.:
            a1 = E
            a2 = np.pi + E
            
            a1diff = np.abs(a_actual - a1)
            a2diff = np.abs(a_actual - a2)
            a_diff_min = min(a1diff, a2diff)
            
            if a_diff_min == a1diff:
                a_calc = a1
            if a_diff_min == a2diff:
                a_calc = a2
                
        return a_calc
        
    def azimuth_alt(x, y, a_actual, a0, x0, y0, A, F, V, S, D, P, Q, e, E):
        r = C * (np.sqrt((x-x0)**2 + (y-y0)**2) + A * (y-y0) * np.cos(F - a0) - A * (x - x0) * np.sin(F - a0))
        u = V * r + S * (np.exp(D * r) - 1.) + P * (np.exp(Q * r * r) - 1.)
        
        b = a0 - E + sc.arctan2((y - y0),(x - x0)) #+ np.pi
       
        z_calc =  np.arccos(np.cos(u) * np.cos(e) - np.sin(u) * np.sin(e) * np.cos(b))
        
        arg = np.sin(b) * np.sin(u) / np.sin(z_calc)
        answ = np.abs(np.arcsin(arg))
        a_calc = determineAz(arg, answ, a_actual, E)
    
        return (a_calc, z_calc, b, u)
    
    
    def minimize(x, y, a, a0, x0, y0, A, F, V, S, D, P, Q, e, E):
        
        for i in range(numSources):
            calc_az_alt[i,0] , calc_az_alt[i,1], w, w = azimuth_alt(x[i], y[i], a[i], a0, x0, y0, A, F, V, S, D, P, Q, e, E) 
        
        F_G = np.subtract(star_az_alt, calc_az_alt)
        F_G = np.abs(F_G) 
        
        F_G_abs = np.abs(F_G)
        finalSum = np.sum(F_G_abs)
        return finalSum
        
    def f():
    #    return a_calc - a_actual
        return np.abs(a_calc - a_actual) * np.sin(z_actual)
        
    def g():
    #    return z_calc - z_actual
        return np.abs(z_calc - z_actual)
        
    def df_da():
        if a_calc > a_actual:
            deriv = np.sin(z_actual)
        if a_calc < a_actual:
            deriv = -np.sin(z_actual)
        else:
            deriv = np.sin(z_actual)
        return deriv
    #
    
    #    return 1.
        
    def dg_dz():
        if z_calc > z_actual:
            deriv = 1.
        if z_calc < z_actual:
            deriv = -1.
        else:
            deriv = 1.
        return deriv
    
    #    return 1.
        
    def da_db():
        return (np.cos(b_calc) * np.sin(u_calc)) / (np.cos(a_calc - E_f) * np.sin(z_calc))
    
    def da_dz():
        return np.negative(np.sin(b_calc) * np.sin(u_calc) * np.cos(z_calc)) / (np.cos(a_calc - E_f) * np.sin(z_calc) * np.sin(z_calc))
    
    def da_du():
        return (np.sin(b_calc) * np.cos(u_calc)) / (np.cos(a_calc - E_f) * np.sin(z_calc))
    
    def da_dE():
        return 1.
        
    def dz_db():
        return np.negative(np.sin(u_calc) * np.sin(e_f) * np.sin(b_calc)) / np.sin(z_calc)
        
    def dz_du():
        return (np.sin(u_calc) * np.cos(e_f) + np.cos(u_calc) * np.sin(e_f) * np.cos(b_calc)) / np.sin(z_calc)
        
    def dz_de():
        return (np.cos(u_calc) * np.sin(e_f) + np.sin(u_calc) * np.cos(e_f) * np.cos(b_calc)) / np.sin(z_calc)
        
    def db_da0():
        return 1.
        
    def db_dx0():
        y_diff = ypix - y0_f
        x_diff_squared = (xpix - x0_f)**2
        return (y_diff / x_diff_squared) * (np.cos(b_calc - a0_f + E_f))**2
        
    def db_dy0():
        x_diff = xpix - x0_f
        return np.negative((np.cos(b_calc - a0_f + E_f))**2) / x_diff
        
    def db_dE():
        return np.negative(1.)
        
    def du_dr():
        l = V_f
        m = S_f * D_f * np.exp(D_f * r)
        n = 2 * Q_f * P_f * r * np.exp(Q_f * r * r)
        return l + m + n
        
    def du_dV():
        return r
        
    def du_dS():
        return np.exp(D_f * r) - 1
    
    def du_dD():
        return S_f * r * np.exp(D_f * r)
        
    def du_dP():
        return np.exp(Q_f * r * r) - 1
        
    def du_dQ():
        return P_f * r * r * np.exp(Q_f * r * r)
        
    def dr_dx0():
        x_diff = xpix - x0_f
        y_diff = ypix - y0_f
        r_simple = np.sqrt(x_diff**2 + y_diff**2)    
        return C * ((-x_diff / r_simple) + A_f * np.sin(F_f - a0_f))
    
    def dr_dy0():
        x_diff = xpix - x0_f
        y_diff = ypix - y0_f
        r_simple = np.sqrt(x_diff**2 + y_diff**2) 
        return C * ((-y_diff / r_simple) - A_f * np.cos(F_f - a0_f))
    
    def dr_da0():
        x_diff = xpix - x0_f
        y_diff = ypix - y0_f
        return A_f * C * (y_diff * np.sin(F_f - a0_f) + x_diff * np.cos(F_f - a0_f))
    
    def dr_dF():    
        x_diff = xpix - x0_f
        y_diff = ypix - y0_f
        return -A_f * C * (y_diff * np.sin(F_f - a0_f) + x_diff * np.cos(F_f - a0_f))
    
    def dr_dA():
        x_diff = xpix - x0_f
        y_diff = ypix - y0_f
        return C * (y_diff * np.cos(F_f - a0_f) - x_diff * np.sin(F_f - a0_f))
            
    def df_da0():
        l = da_db() * db_da0()
        m = da_du() * du_dr() * dr_da0()
        n = da_dz() * (dz_du() * du_dr() * dr_da0() + dz_db() * db_da0())
        return df_da() * (l + m + n)
        
    def dg_da0():
        a = dz_du() * du_dr() * dr_da0()
        b = dz_db() * db_da0()
        return dg_dz() * (a + b)
        
    def df_dx0():
        l = da_db() * db_dx0()
        m = da_du() * du_dr() * dr_dx0()
        n = da_dz() * (dz_du() * du_dr() * dr_dx0() + dz_db() * db_dx0())
        return df_da() * (l + m + n)
        
    def df_dy0():
        l = da_db() * db_dy0()
        m = da_du() * du_dr() * dr_dy0()
        n = da_dz() * (dz_du() * du_dr() * dr_dy0() + dz_db() * db_dy0())
        return df_da() * (l + m + n)
        
    def dg_dx0():
        l = dz_du() * du_dr() * dr_dx0()
        m = dz_db() * db_dx0()
        return dg_dz() * (l + m)
        
    def dg_dy0():
        l = dz_du() * du_dr() * dr_dy0()
        m = dz_db() * db_dy0()
        return dg_dz() * (l + m)
    
    def df_dV():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m)
        return base * du_dV()
        
    def dg_dV():
        base = dg_dz() * dz_du()
        return base * du_dV()
        
    def df_dS():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m)
        return base * du_dS()
        
    def dg_dS():
        base = dg_dz() * dz_du()
        return base * du_dS()
    
    def df_dD():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m)
        return base * du_dD()
        
    def dg_dD():
        base = dg_dz() * dz_du() * du_dr()
        return base * du_dD()
        
    def df_dP():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m)
        return base * du_dP()
        
    def dg_dP():
        base = dg_dz() * dz_du()
        return base * du_dP()
    
    def df_dQ():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m)
        return base * du_dQ()
        
    def dg_dQ():
        base = dg_dz() * dz_du()
        return base * du_dQ()
        
    def df_de():
        return df_da() * da_dz() * dz_de()
        
    def dg_de():
        return dg_dz() * dz_de()
        
    def df_dE():
        l = da_dE()
        m = (da_db() + da_dz() * dz_db()) * db_dE()
        return df_da() * (l + m)
        
    def dg_dE():
        return dg_dz() * dz_db() * db_dE()
        
    def df_dA():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m) * du_dr()
        return base * dr_dA()
        
    def dg_dA():
        base = dg_dz() * dz_du() * du_dr()
        return base * dr_dA()
        
    def df_dF():
        l = da_du()
        m = da_dz() * dz_du()
        base = df_da() * (l + m) * du_dr()
        return base * dr_dF()
        
    def dg_dF():
        base = dg_dz() * dz_du() * du_dr()
        return base * dr_dF()
        
    
    a = star_az_alt[:,0]    
    z = star_az_alt[:,1]
    x = star_x_y[:,0]
    y = star_x_y[:,1]
    calc_az_alt = np.zeros(shape=star_az_alt.shape)
#    print star_az_alt
    
    
    #initialize parameters
    a0 = -1.           #angle between the x axis and the direction to the south
    x0 = 315.             #center of x-coord projection on xy plane
    y0 = 250.               #center of y-coord projection on xy plane
    A = 0.0                 #related to eccentricity of proj due to tilt of camera
    F = -15                 #related to direction of inclination of projection due to tilt of camera
    C = 1.                  #global scale factor (shift of plate along the axis)
    V = 0.72 * mm2pix       #lens constant
    S = 0.0009              #lens constant
    D = 0.6 * mm2pix        #lens constant
    P = 0.0                 #lens constant
    Q = 0.0 * mm2pix        #lens constant
    e = .002
    
#    a0 = 1.51174644e+01
#    x0 = 3.37437442e+02
#    y0 = 2.08128479e+02
#    A = 4.87413308e-02
#    F = -3.15101264e+00
#    V = 5.33041484e-03
#    S = 2.38514600e-03
#    D = -1.06291427e+00
#    P = -5.98757680e+04
#    Q = 1.31505071e-11
#    e = 1.38395489e-01
#    E = 1.64685925e-01
#    return [a0, x0, y0, A, F, V, S, D, P, Q]
#    return [a0, x0, y0, A, F, V, S, D, P, Q]

    #e_E = .1
       
    #a0_i = 2.21550496e+00
    #x0_i =  3.46053259e+02
    #y0_i = 2.16005659e+02
    #A_i =-5.23674642e-02
    #F_i = -1.90523719e+01
    #V_i =  5.34564287e-03
    #S_i = 5.04970555e-03
    #D_i =  -3.30385847e+01
    #P_i = -7.95110819e+05
    #Q_i =  9.05285346e-13
    #e =6.67688137e-02
    #E_mi =  3.21803481e-03

    #Fit 1
    initialParams_az = (a0, x0, y0)
    
    paramFits_az, cov_az = op.curve_fit(az_initial_fit, star_x_y, a, p0=initialParams_az)    
    
    print 'Parameter Fits, Round 1: '
    print paramFits_az
    print ''
    
    #Fit 2
    initialParams_alt = (paramFits_az[0], paramFits_az[1], paramFits_az[2], A, F, V, S, D, P, Q)
    
    paramFits_alt, cov_alt = op.curve_fit(alt_initial_fit, star_x_y, z, p0=initialParams_alt)
        
    print 'Parameter Fits, Round 2: '
    print paramFits_alt
    print ''
          
    #Name initial params
    a0_i = paramFits_alt[0]
    x0_i = paramFits_alt[1]
    y0_i = paramFits_alt[2]
    A_i = paramFits_alt[3]
    F_i = paramFits_alt[4]
    V_i = paramFits_alt[5]
    S_i = paramFits_alt[6]
    D_i = paramFits_alt[7]
    P_i = paramFits_alt[8]
    Q_i = paramFits_alt[9]
#    return [a0_i, x0_i, y0_i, A_i, F_i, V_i, S_i, D_i, P_i, Q_i]

    #Find best initial guess for E    
    E_tries = np.linspace(-4*np.pi, 4*np.pi, 500)
    E_tries_min = np.zeros(shape=E_tries.shape)
    numTries = E_tries.shape[0]
    
    
    for i in range(numTries):
        x_new = x0_i + (e / V_i) * np.cos(E_tries[i])
        y_new = y0_i + (e / V_i) * np.sin(E_tries[i])
        E_tries_min[i] = minimize(x, y, a, a0_i, x_new, y_new, A_i, F_i, V_i, S_i, D_i, P_i, Q_i, e, E_tries[i])    
       
    E_min_index = np.argmin(E_tries_min)
    E_min = E_tries[E_min_index] #a0_i=-1.28440609e+00
    
    #Recalculate x0 and y0 guesses based on new E_min
    x0_i = x0_i + (e / V_i) * np.cos(E_min)
    y0_i = y0_i + (e / V_i) * np.sin(E_min)
    return [a0_i, x0_i, y0_i, A_i, F_i, V_i, S_i, D_i, P_i, Q_i, e, E_min]

    
    #a0_i =  2.21534040e+00
    #x0_i = 3.46053910e+02
    #y0_i = 2.16005157e+02
    #A_i =-5.23690633e-02
    #F_i = -1.90525373e+01
    #V_i =  5.34564874e-03
    #S_i =  5.05019467e-03
    #D_i =  -3.30385847e+01
    #P_i = -7.95110819e+05
    #Q_i = 9.05305699e-13
    #e =6.67687719e-02
    #E_min =  3.13580778e-03
    #
    #a0_i = 2.21550496e+00
    #x0_i =  3.46053259e+02
    #y0_i = 2.16005659e+02
    #A_i =-5.23674642e-02
    #F_i = -1.90523719e+01
    #V_i =  5.34564287e-03
    #S_i = 5.04970555e-03
    #D_i =  -3.30385847e+01
    #P_i = -7.95110819e+05
    #Q_i =  9.05285346e-13
    #e =6.67688137e-02
    #E_min =  3.21803481e-03
    #
    #a0_i = 2.21550496e+00
    #x0_i =  3.46053259e+02
    #y0_i = 2.16005659e+02
    #A_i =-5.23674642e-02
    #F_i = -1.90523719e+01
    #V_i =  5.34564287e-03
    #S_i = 5.04970555e-03
    #D_i =  -3.30385847e+01
    #P_i = -7.95110819e+05
    #Q_i =  9.05285346e-13
    #e =6.67688137e-02
    #E_min =  3.21803481e-03
    a0_i = 1.51174644e+01
    x0_i = 3.37437442e+02
    y0_i = 2.08128479e+02
    A_i = 4.87413308e-02
    F_i = -3.15101264e+00
    V_i = 5.33041484e-03
    S_i = 2.38514600e-03
    D_i = -1.06291427e+00
    P_i = -5.98757680e+04
    Q_i = 1.31505071e-11
    e = 1.38395489e-01
    E_min = 1.64685925e-01

    #Calculate az/alt with set of parameters
    for i in range(numSources):
            calc_az_alt[i,0] , calc_az_alt[i,1], w, w = azimuth_alt(x[i], y[i], a[i], a0_i, x0_i, y0_i, A_i, F_i, V_i, S_i, D_i, P_i, Q_i, e, E_min) 
    print ''
    print ''
           
    print '********SET OF INITIAL PARAMETERS*********'    
    print 'a0: %0.4f' % a0_i
    print 'x0: %0.4f' % x0_i
    print 'y0: %0.4f' % y0_i
    print 'A: %0.4f' % A_i
    print 'F: %0.4f' % F_i
    print 'V: %0.4f' % V_i
    print 'S: %0.4f' % S_i
    print 'D: %0.4f' % D_i
    print 'P: %0.4f' % P_i
    print 'Q: %0.4f' % Q_i
    print 'e: %0.4f' % e
    print 'E: %0.4f' % E_min
    print ''
    
    
    #Plot alt
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.scatter(x, y, z, c='r')
    
    ax1.scatter(x, y, calc_az_alt[:,1], c = 'b')
    
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Altitude')
    
    
    #Plot az
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(x, y, a, c='r')
    
    ax2.scatter(x, y, calc_az_alt[:,0], c = 'b')
    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Azimuth')
    
    plt.show()  
    
    
    
    #parameter tolerances
    da0 = 0.00001
    dx0 = 0.00001
    dy0 = 0.00001
    dA  = 0.00001
    dF  = 0.00001
    dV  = 0.00001
    dS  = 0.00001
    dD  = 0.00001
    dP  = 0.00001
    dQ  = 0.00001
    de  = 0.00001
    dE  = 0.00001
    
    
    
    Params = np.array([[a0_i], [x0_i], [y0_i], [A_i], [F_i], [V_i], [S_i], [D_i], [P_i], [Q_i], [e], [E_min]])
    
    numParams = Params.shape[0]
    
    derivativeMatrix = np.zeros(shape=(2*numSources, numParams))
    
    dp = np.zeros(shape = (numParams, 1))
    
    tolerance_dp = np.negative(np.ones(shape = (numParams, 1)))
    
    f_g_array = np.zeros(shape = (2*numSources, 1))
    
    
    while np.any(np.abs(dp) > tolerance_dp):
    
        
        tolerance_dp[0], tolerance_dp[1], tolerance_dp[2], tolerance_dp[3], tolerance_dp[4], tolerance_dp[5], tolerance_dp[6], tolerance_dp[7], tolerance_dp[8], tolerance_dp[9], tolerance_dp[10], tolerance_dp[11] = da0, dx0, dy0, dA, dF, dV, dS, dD, dP, dQ, de, dE
        
        print ''
        print 'NEW ROUND'
        print ''
        print 'DP'    
        print dp
            
        for i in range(numParams):
            Params[i] += 0.0001*dp[i]
                
        a0_f, x0_f, y0_f, A_f, F_f, V_f, S_f, D_f, P_f, Q_f, e_f, E_f = Params
        
        print''
        print 'PARAMS'
        print Params
        
        for i in range(numSources):
            
            xpix = star_x_y[i,0]
            ypix = star_x_y[i,1]
            a_actual = star_az_alt[i,0]
            z_actual = star_az_alt[i,1]
            a_calc, z_calc, b_calc, u_calc = azimuth_alt(xpix, ypix, a_actual, a0_f, x0_f, y0_f, A_f, F_f, V_f, S_f, D_f, P_f, Q_f, e_f, E_f)
            r = C * (np.sqrt((xpix - x0_f)**2 + (ypix - y0_f)**2) + A_f * (ypix - y0_f) * np.cos(F_f - a0_f) - A_f * (xpix-x0_f) * np.sin(F_f - a0_f))
            
            derivativeMatrix[i,0] = df_da0()
            derivativeMatrix[i,1] = df_dx0()
            derivativeMatrix[i,2] = df_dy0()
            derivativeMatrix[i,3] = df_dA()
            derivativeMatrix[i,4] = df_dF()
            derivativeMatrix[i,5] = df_dV()
            derivativeMatrix[i,6] = df_dS()
            derivativeMatrix[i,7] = df_dD()
            derivativeMatrix[i,8] = df_dP()
            derivativeMatrix[i,9] = df_dQ()
            derivativeMatrix[i,10] = df_de()
            derivativeMatrix[i,11] = df_dE()
            
            f_g_array[i] = -f()  
            
        for i in range(numSources, 2*numSources):
        
            xpix = star_x_y[i-numSources,0]
            ypix = star_x_y[i-numSources,1]
            a_actual = star_az_alt[i-numSources,0]
            z_actual = star_az_alt[i-numSources,1]
            a_calc, z_calc, b_calc, u_calc = azimuth_alt(xpix, ypix, a_actual, a0_f, x0_f, y0_f, A_f, F_f, V_f, S_f, D_f, P_f, Q_f, e_f, E_f)
            r = C * (np.sqrt((xpix - x0_f)**2 + (ypix - y0_f)**2) + A_f * (ypix - y0_f) * np.cos(F_f - a0_f) - A_f * (xpix - x0_f) * np.sin(F_f - a0_f))
            
            
            derivativeMatrix[i,0] = dg_da0()
            derivativeMatrix[i,1] = dg_dx0()
            derivativeMatrix[i,2] = dg_dy0()
            derivativeMatrix[i,3] = dg_dA()
            derivativeMatrix[i,4] = dg_dF()
            derivativeMatrix[i,5] = dg_dV()
            derivativeMatrix[i,6] = dg_dS()
            derivativeMatrix[i,7] = dg_dD()
            derivativeMatrix[i,8] = dg_dP()
            derivativeMatrix[i,9] = dg_dQ()
            derivativeMatrix[i,10] = dg_de()
            derivativeMatrix[i,11] = dg_dE()
            
            f_g_array[i] = -g()
        
        dp, residuals, rank, s = np.linalg.lstsq(derivativeMatrix, f_g_array)
    
    
    #Calculate az/alt with set of final parameters
    for i in range(numSources):
        calc_az_alt[i,0] , calc_az_alt[i,1], w, w = azimuth_alt(x[i], y[i], a[i], Params[0], Params[1], Params[2], Params[3], Params[4], Params[5], Params[6], Params[7], Params[8], Params[9], Params[10], Params[11]) 
    
    
       
    ##Plot alt
    #fig1 = plt.figure()
    #ax1 = fig1.add_subplot(111, projection='3d')
    #ax1.scatter(x,y,z, c='r')
    #
    #ax1.scatter(x,y,calc_az_alt[:,1], c = 'b')
    #
    #ax1.set_xlabel('X')
    #ax1.set_ylabel('Y')
    #ax1.set_zlabel('Altitude')
    #
    #
    ##Plot az
    #fig2 = plt.figure()
    #ax2 = fig2.add_subplot(111, projection='3d')
    #ax2.scatter(x,y,a, c='r')
    #
    #ax2.scatter(x,y,calc_az_alt[:,0], c = 'b')
    #
    #ax2.set_xlabel('X')
    #ax2.set_ylabel('Y')
    #ax2.set_zlabel('Azimuth')
    #
    #plt.show()  
    #
    #Error in final az/alt values
    diff_az_alt = (star_az_alt - calc_az_alt)
    
    print ''
    
    
    print '********FINAL PARAMETERS AND ERRORS*********'    
    print 'a0:  %3.6f, %.2e' % (Params[0], dp[0])
    print 'x0:  %3.6f, %.2e' % (Params[1], dp[1])
    print 'y0:  %3.6f, %.2e' % (Params[2], dp[2])
    print 'A:   %3.6f, %.2e' % (Params[3], dp[3])
    print 'F:   %3.6f, %.2e' % (Params[4], dp[4])
    print 'V:   %3.6f, %.2e' % (Params[5], dp[5])
    print 'S:   %3.6f, %.2e' % (Params[6], dp[6])
    print 'D:   %3.6f, %.2e' % (Params[7], dp[7])
    print 'P:   %3.6f, %.2e' % (Params[8], dp[8])
    print 'Q:   %3.6f, %.2e' % (Params[9], dp[9])
    print 'e:   %3.6f, %.2e' % (Params[10], dp[10])
    print 'E:   %3.6f, %.2e' % (Params[11], dp[11])
    print ''
    
    
    print ''
    print ''
    print '*******RESIDUALS OF AZ ALT*******'
    print diff_az_alt
    
    #plt.plot(x, diff_az_alt[:,0], 'r+')
    #plt.plot(y, diff_az_alt[:,0], 'b+')
    #plt.xlabel('x (red) and y (blue)')
    #plt.ylabel('Azimuth')
    #plt.show()
    #
    #plt.plot(x, diff_az_alt[:,1], 'r+')
    #plt.plot(y, diff_az_alt[:,1], 'b+')
    #plt.xlabel('x (red) and y (blue)')
    #plt.ylabel('Altitude')
    #plt.show()
    return Params
    
if __name__ == '__main__':
    main()