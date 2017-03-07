'''
CluStR takes a FITS cluster catalog and calculates various scaling relations
using either the Kelly method (using `linmix`) or the Mantz method (using `lrgs`).
Modifies `devon_scaling_relations` by Devon Hollowood.

Spencer Everett, UCSC, 3/2017

Initially based upon `devon_scaling_relations` by Devon Hollowood.
'''

#pylint: disable=invalid-name
#pylint: disable=no-member

import argparse
import os, corner, ast
import astropy.io.fits as fits
import numpy as np
import reglib # Regression library
import plotlib # Plotting library
import PyPDF2

# Parameter list used throughout. See `param.config`
parameters = {}

def Ez(z,Om=0.3,H_0=0.7):
    ''' Calculate E(z) for Om=0.3, H_0=0.7 cosmology. '''
    h = H_0/100.
    return np.sqrt(Om*(1.+z)**3 + h)

def fits_label(axis_name):
    ''' Get the FITS column label for `axis_name` '''
    labels = {
        'lambda':'lambda',
        'l500kpc':'500_kiloparsecs_band_lumin',
        'lr2500':'r2500_band_lumin',
        'lr500':'r500_band_lumin',
        'lr500cc':'r500_core_cropped_band_lumin',
        't500kpc':'500_kiloparsecs_temperature',
        'tr2500':'r2500_temperature',
        'tr500':'r500_temperature',
        'tr500cc':'r500_core_cropped_temperature'
    }
    return labels[axis_name]

def check_flag(flag):
    ''' Checks whether the flag type is boolean, cutoff, or a range. '''
    #FIX: Finish!!

    boolean = {
        'Analyzed',
        'Detected',
        'merger',
        'bad_redmapper_pos',
        'bad_xray_pos',
        'bad_pos_other',
        'on_chip_edge',
        'edge_exclude_centering',
        'edge_exclude_r2500',
        'edge_exclude_r500',
        'off_axis_chip',
        'serendipitous',
        'overlap_r2500',
        'overlap_r500'
    }

    cutoff = {
        'offset_r500',
        'offset_r2500'
    }

    # Mispelled as 'range' is reserved
    rnge = {
        'Redshift',
        'redMaPPer_ra',
        'redMaPPer_dec',
        'r500_ra',
        'r500_dec',
        'r2500_ra',
        'r2500_dec',
        '500_kiloparsecs_ra',
        '500_kiloparsecs_dec',
        'lambda',
        'r500_band_lumin',
        'r500_temperature',
        'r500_core_cropped_band_lumin',
        'r500_core_cropped_temperature',
        'r2500_band_lumin',
        'r2500_temperature',
        '500_kiloparsecs_band_lumin',
        '500_kiloparsecs_temperature'
    }

    if flag in boolean:
        return 'bool'
    elif flag in cutoff:
        return 'cutoff'
    elif flag in rnge:
        return 'range'
    else:
        print('WARNING: {} is an invalid flag entry. Ignoring flag.'.format(flag))

    return None

def get_data(options):
    ''' Get x, y, x errors, y errors '''

    hdulist = fits.open(options.catalog)
    data = hdulist[1].data

    label_x = fits_label(options.x)
    label_y = fits_label(options.y)
    x = data[label_x]
    y = data[label_y]

    # Numer of original data
    N = np.size(x)

    # divide luminosities by e[z]
    # FIX: Probably shouldn't hardcode this!
    if label_x[0] == 'l' and label_y != 'lambda':
        x /= Ez(data['redshift'])
    if label_y[0] == 'l' and label_x != 'lambda':
        y /= Ez(data['redshift'])

    x_err = (data[label_x+'_err_low'] + data[label_x+'_err_high']) / 2.
    y_err = (data[label_y+'_err_low'] + data[label_y+'_err_high']) / 2.

    # Now take out any flagged data
    flags = options.flags
    if flags is not None:
        # FIX: Should be more error handling than this!
        # FIX: Should write a method to ensure all the counts are what we expect
        for flag in flags:
            flag_type = check_flag(flag)

            if flag_type == 'bool':
                # If flag is boolean, remove those that are True
                flagged = np.where(data[flag] == True)
                x[flagged] = -1
                y[flagged] = -1
                L = np.size(flagged)
                print('Removed {} data points due to {} flag'.format(L,flag))
            elif flag_type == 'cutoff':
                # Remove data above/below desired cutoff
                pass #FIX
            elif flag_type == 'range':
                # Remove data outside of inputted range
                pass # FIX
            else:
                print('WARNING: Inputted flag {} is not valid. Ignoring flag for analysis.'.format(flag))
                continue

    # Take rows with good data, and all flagged data removed
    good_rows = np.all([x != -1, y != -1], axis=0)

    x = x[good_rows]
    y = y[good_rows]
    x_err = x_err[good_rows]
    y_err = y_err[good_rows]

    print('Accepted {} data out of {}'.format(np.size(x),N))

    print 'mean x error:', np.mean(x_err)
    print 'mean y error:', np.mean(y_err)

    hdulist.close()

    # return shifted x, x pivot, y, x errors, y errors
    return (x, y, x_err, y_err)

def scale(x, y, x_err, y_err):
    ''' Scale data for fitting '''
    log_x = np.log(x)
    x_piv = np.median(log_x)
    log_y = np.log(y)
    return (log_x-x_piv, log_y, x_err/x, y_err/y, x_piv)

def unscale(x, y, x_err, y_err, x_piv):
    ''' Recover original data from fit-scaled data '''
    return (np.exp(x + x_piv), np.exp(y), x_err * x, y_err * y)

def scaled_fit_to_data(x_min, x_max, x_piv, scaled_fit):
    ''' Get a data set from a scaled fit '''
    (fit_int, fit_slope, fit_sig) = scaled_fit
    scaled_x = np.linspace(x_min, x_max, 101)
    scaled_y = np.mean(fit_int) + np.mean(fit_slope) * scaled_x
    scaled_x_errs = np.zeros(101)
    scaled_y_errs = np.ones(101)*np.mean(fit_sig)
    unscaled_data = unscale(scaled_x, scaled_y, scaled_x_errs, scaled_y_errs, x_piv)
    return unscaled_data

def fit(method,x_obs, y_obs, x_err, y_err,nmc=5000):
    ''' Calculates fit using the Kelly (linmix) and/or Mantz (lrgs) methods and returns
        their respective markov chains for the intercept, slope, and sigma. Does only
        Kelly method by default.'''

    #print "median log(x-x_piv) =", np.median(x_obs)
    assert np.median(x_obs) == 0.0

    methods = return_methods_list(method)

    # Set default parameter markov chains to None in case both methods aren't used
    kelly_b, kelly_m, kelly_sig = None, None, None
    mantz_b, mantz_m, mantz_sig = None, None, None

    # Iterate through desired methods
    for method in methods:
        if method == 'kelly':
            print "Using Kelly Algorithm..."
            kelly_b, kelly_m, kelly_sig = reglib.run_linmix(x_obs, y_obs, x_err, y_err)

        if method == 'mantz':
            print "Using Mantz Algorithm..."
            mantz_b, mantz_m, mantz_sig = reglib.run_lrgs(x_obs,y_obs,x_err,y_err, nmc=nmc)

    return (kelly_b, kelly_m, kelly_sig), (mantz_b, mantz_m, mantz_sig)

def return_methods_list(method):
    ''' Used to convert a method input choice into a list of used method types.'''
    if method is None:
        # Default is Kelly method
        methods = ['kelly']
    elif method.lower() == 'both':
        methods = ['kelly', 'mantz']
    elif method.lower() == 'kelly' or method.lower() == 'mantz':
        # Needed for loop structure
        methods = [method]
    else:
        print "WARNNG: Only `kelly`, `mantz`, or `both` are valid method options. Will use Kelly method instead."
        methods = ['kelly']

    return methods

def set_parameters(file):
    ''' Set useful parameters from config file'''

    config_file = open(file)
    for line in config_file:
        # Ignore empty lines and comments:
        if line[0:2] == '\n': continue
        if line[0] == '#': continue
        line.strip()
        line = line.split('#')[0]
        # Remove whitespace and interpret Name:Value pairs:
        line = ''.join(line.split())
        line = line.split(':')
        Name, Value = line[0], line[1]
        parameters[Name] = Value

    config_file.close()

    return

def parse_opts():
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser()
    # Required argument for catalog
    parser.add_argument('catalog', help='FITS catalog to open')
    # Required arguement for axes
    valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500', 'tr500', 'tr500cc',
                  'lambda']
    parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
    parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
    # Optional argument for file prefix
    parser.add_argument('-p', '--prefix', help='prefix for output file')
    # Optional argument for regression method(s)
    methods = ['kelly','mantz','both']
    parser.add_argument('-m', '--method',help='Choose the `kelly` or `mantz` regression method (or `both`)', choices=methods)
    # Optional arguments for any flag cuts
    # FIX: in the future, make an allowed choices vector work!
    parser.add_argument('-f', '--flags',nargs='+', type=str,help='Input any desired flag cuts as a list of flag names (with "" and no spaces!)')
    # Optional argument for which files to be saved
    # FIX: Implement!

    return parser.parse_args()

def check_dependencies():
    '''
    In the future, this function will check if all required packages are installed, and, if not, ask
    if these packages should be downloaded. Otherwise exists program.
    '''
    # FIX: Implement!
    # linmix
    # lrgs
    # corner
    # pypdf2
    # check for others!
    pass

def main(): #pylint: disable=missing-docstring
    print('Checking dependencies...')
    check_dependencies()

    # Parse all inputted options, regardless of param.config file
    options = parse_opts()

    # Set useful parameters from configure file
    config_file = 'param.config'
    set_parameters(config_file)

    # Determine methods list
    method = options.method

    # Determine flags
    flags = options.flags

    print('Inputted options:',options)
    print('Grabbing data...')

    # Grab and process data from catalog, including flag removal
    data_obs = get_data(options)

    # Scale for linear fitting
    scaled_data = scale(*data_obs)

    print('Fitting data...')

    # Fit data using linmix, lrgs, or both
    kelly_scaled_fit, mantz_scaled_fit = fit(method,*scaled_data[:4])
    (x_min, x_max) = (np.min(scaled_data[0]), np.max(scaled_data[0]))

    print('Making plots...')

    # Make all desired plots
    plotlib.make_plots(options, data_obs, kelly_scaled_fit, mantz_scaled_fit, scaled_data[4], x_min, x_max)

    print('Done!')

if __name__ == '__main__':
    main()
