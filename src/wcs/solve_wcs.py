'''
solve_wcs.py - Functions to solve image x, y to ra, dec.

Author: 
    Elisabeth Heinrich-Josties (eheinrich-josties@lcogt.net)

May 2015

'''
import logging
reload(logging)

import calculation_steps as steps
from utils.astro_utils import dist_to_pix
from derivatives import Derivatives
from wcs.solved_parameters import parameters

logger = logging.getLogger(__name__)



def solve(data_filename, image_filename, ra_col, dec_col, x_col, y_col, fraction=0.001,
          param_tolerance=0.00001, initial_params=parameters('generic'), redo_initial_fits=True):
    '''
    Solves for the projection parameters of a given camera, which will allow accurate 
    world coordinates to be determined from x,y pixel position. The projection equations
    and parameters are described by Borovicka et al. and can be found 
    here: http://adsabs.harvard.edu/abs/1995A%26AS..112..173B
    
    Note: The solving process is very sensitive to the initial parameters supplied.
    A reasonable set of input parameters is used to begin with, but it is likely that
    they will have to be tweak before all of the equations can solve properly.
    
    I suggest running the program as follows:
    - run solve with only the required arguments. If you get something like:
                ValueError: Unexpected value for arg, got [ nan].
      Then reduce fraction to something really small, like 0.0000000001,
      and let that run until it converges. It might take a while. Just kep an eye 
      on the numbers that it outputs during the loops (dp array) to make sure they are steadily
      incrementing--the values should not be jumping around. 
    - After getting to a solution, run solve again on the output final_params, with 
      same or larger fraction, and lower tolerances. If it seems like the parameter errors
      arent changing, try a larger fraction.
    - Once you have sufficiently small dp/residuals, and the outputted graphs look reasonable,
      the solved_parameters module can be modified to include the new parameter values.
      
    -OR, email me!
      
    :param data_filename: Path to the text file that contains data of x,y pixel positions and
                          their corresponding ra,dec. This file should be made by hand, and 
                          there should be at least 20 data points.
    :param image_filename: Path to the FITS file from which the x,y and ra,dec pairs are found.
    :param ra_col: Column number of text file that has RAs, zero-indexed
    :param dec_col: Column number of text file that has decs, zero-indexed
    :param x_col: Column number of text file that has x pixel positions, zero-indexed
    :param y_col: Column number of text file that has y pixel positions, zero-indexed
    :param fraction: float multiplier to determine step size in the linearization process
    :param param_tolerance: float, the smaller, the smaller the errors on the parameters must be before convergence
    :param initial_params: Can supply a dictionary of initial projection parameters for the code to start with
    :param redo_initial_fits: Boolean, can set it to false to skip the initial fits. This is useful for if parameters
                              that are close are already known and supplied to the code--setting this to false will
                              ensure that those parameters arent modified before they get to the linearization step
    :return: Final fitted parameters as a dictionary, and residuals (calculated az/alt - measured az/alt in pixels)
    
    If the code has been looping for a while, and you wish to see how the results are looking, Ctrl+C, and the current 
    parameters will be returns, and plots of real and calculated azimuth/ zenith distance will be outputted. You can then set
    initial_params=final_params, and redo_initial_fits=False if the results look like theyre starting to converge.
    
    '''
    logging.basicConfig(level=logging.INFO, format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')    
    
    logger.info('Gathering data from %s' % data_filename)
    
    star_ra_dec, star_x_y = steps.get_data(data_filename, ra_col, dec_col, x_col, y_col)

    star_az_zendist = steps.az_zendist_of_data(image_filename, star_ra_dec)
    
    if redo_initial_fits:
        # Do initial fits of u and b equations to get better initial guess on 
        # all parameters except for E and e. Must reassign initial params 
        initial_params = steps.initial_fits(star_x_y, star_az_zendist, initial_params)
    
        # Find better initial guess for E, x0, y0 and reassign initial_params
        initial_params = steps.guess_for_E_x0_y0(star_x_y, star_az_zendist, initial_params)
    
    # linearize and solve for final parameters
    derivatives_to_solve = Derivatives(initial_params)
    
    try:
        final_params, dp = derivatives_to_solve.leastSqSolve(star_x_y, star_az_zendist, fraction=fraction, param_tolerance=param_tolerance)
    except KeyboardInterrupt:
        calc_az_zendist, residuals = calc_az_zendist, residuals = steps.calc_az_zendist_with_params(star_x_y, star_az_zendist, final_params)
        steps.interactive_3D_plot(star_x_y, star_az_zendist, calc_az_zendist)
        return final_params, dp, residuals
        
    # Calculate az/alt with set of final parameters and plot
    calc_az_zendist, residuals = steps.calc_az_zendist_with_params(star_x_y, star_az_zendist, final_params)
    steps.interactive_3D_plot(star_x_y, star_az_zendist, calc_az_zendist)

    return final_params, dp, dist_to_pix(residuals, pixscale=0.3, u='radian')
    
    
    