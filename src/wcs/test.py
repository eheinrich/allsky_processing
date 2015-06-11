"""
Tests to make sure refactoring is properly done.

@author: eheinrich
"""
import sys
sys.path.append('/home/eheinrich/Desktop/Code/allsky_processing/src')
import wcs_paper
import solve_wcs
import pytest

def test():
    image_data = '/home/eheinrich/Desktop/Code/CloudsCode/china/Inputs/stars_china.dat'
    current_image = '/home/eheinrich/Desktop/Code/CloudsCode/china/Inputs/allsky_china.fits'
    
    original_method_params = wcs_paper.main()
    
    new_method_params = solve_wcs.solve(image_data, current_image, 1, 2, 3, 4)  
    print original_method_params
    print ''
    print new_method_params
    
#    return original_method_params, new_method_params

    assert original_method_params[0] == new_method_params['a0'] #, 'a0 diff in param: %f' % original_method_params[0]-new_method_params['a0']
    assert original_method_params[1] == new_method_params['x0'] #, 'x0 diff in param: %f' % original_method_params[0]-new_method_params['x0']
    assert original_method_params[2] == new_method_params['y0'] #, 'y0 diff in param: %f' % original_method_params[0]-new_method_params['y0']
    assert original_method_params[3] == new_method_params['A'] #, 'A diff in param: %f' % original_method_params[0]-new_method_params['A']
    assert original_method_params[4] == new_method_params['F'] #, 'F diff in param: %f' % original_method_params[0]-new_method_params['F']
    assert original_method_params[5] == new_method_params['V'] #, 'V diff in param: %f' % original_method_params[0]-new_method_params['V']
    assert original_method_params[6] == new_method_params['S'] #, 'S diff in param: %f' % original_method_params[0]-new_method_params['S']
    assert original_method_params[7] == new_method_params['D'] #, 'D diff in param: %f' % original_method_params[0]-new_method_params['D']
    assert original_method_params[8] == new_method_params['P'] #, 'P diff in param: %f' % original_method_params[0]-new_method_params['P']
    assert original_method_params[9] == new_method_params['Q'] #, 'Q diff in param: %f' % original_method_params[0]-new_method_params['Q']
    assert original_method_params[10] == new_method_params['e'] #, 'e diff in param: %f' % original_method_params[0]-new_method_params['e']
    assert original_method_params[11] == new_method_params['E'] #, 'E diff in param: %f' % original_method_params[0]-new_method_params['E']

    return original_method_params, new_method_params
    
if __name__ == '__main__':
    test()