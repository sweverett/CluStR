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
import astropy.io.fits as fits
import numpy as np
import reglib # Regression library
import plotlib # Plotting library

# Parameter list used throughout. See `param.config`
PARAMETERS = {}

# Method List
METHODS = []

def Ez(z):
    ''' Calculate E(z) for Om=0.3, H_0=0.7 cosmology. '''
    Om = PARAMETERS['Om']
    H_0 = PARAMETERS['H_0']
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
    #FIX: Flags are now more versatile, so change up the cutoff and

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

    # List of possible cutoff and range flags
    cut_or_range = {
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
        '500_kiloparsecs_temperature',
        'offset_r500',
        'offset_r2500'
    }

    if flag in boolean:
        bool_type = PARAMETERS[flag+'_bool_type']
        if bool_type is True or bool_type is False:
            return 'bool'
        else:
            print('WARNING: {} is an invalid type entry. Ignoring {} flag.'.format(cut_type,flag))
            return None

    elif flag in cut_or_range:
        # Tries cutoff, goes to range if fails. Can only pick one in config file!
        try:
            cut_type = PARAMETERS[flag+'_cut_type'].lower()
            if cut_type == 'above' or cut_type == 'below':
                return 'cutoff'
            else:
                print('WARNING: {} is an invalid type entry. Ignoring {} flag.'.format(cut_type,flag))
                return None

        except:
            try:
                range_min = PARAMETERS[flag+'_range_min']
                range_max = PARAMETERS[flag+'_range_max']
                range_type = PARAMETERS[flag+'_range_type']
                if range_type == 'inside' or range_type = 'outside':
                    return 'range'
                else:
                    print('WARNING: {} is an invalid type entry. Ignoring {} flag.'.format(range_type,flag))

            except:
                # FIX!!!
                pass
                
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
                # Remove data that is flagged True or False
                bool_type = PARAMETERS[flag+'_bool_type']
                if bool_type is True or bool_type is False:
                    # Want to remove the complement of bool_type
                    flagged = np.where(data[flag] == (not bool_type) )
                    x[flagged] = -1
                    y[flagged] = -1
                    L = np.size(flagged)
                    print('Removed {} data points due to {} flag of type {}'.format(L,flag,flag_type))
                else:
                    print('WARNING: Boolean type must be `true` or `false`. Ignoring {} flag.'.format(flag))
                    continue
            elif flag_type == 'cutoff':
                # Remove data above/below desired cutoff
                cutoff = PARAMETERS[flag+'_cut']
                if PARAMETERS[flag+'_cut_type'].lower() == 'above':
                    flagged = np.where(data[flag] > cutoff)
                elif PARAMETERS[flag+'_cut_type'].lowe() == 'below:':
                    flagged = np.where(data[flag] < cutoff)
                else:
                    print('WARNING: Cutoff type must be `above` or `below`. Ignoring {} flag.'.format(flag))
                    continue
                x[flagged] = -1
                y[flagged] = -1
                L = np.size(flagged)
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

def fit(x_obs, y_obs, x_err, y_err,nmc=5000):
    ''' Calculates fit using the Kelly (linmix) and/or Mantz (lrgs) methods and returns
        their respective markov chains for the intercept, slope, and sigma. Does only
        Kelly method by default.'''

    print "median log(x-x_piv) =", np.median(x_obs)
    #assert np.median(x_obs) == 0.0

    # Set default parameter markov chains to None in case both methods aren't used
    kelly_b, kelly_m, kelly_sig = None, None, None
    mantz_b, mantz_m, mantz_sig = None, None, None

    # Iterate through desired methods
    for method in METHODS:
        if method == 'kelly':
            print "Using Kelly Algorithm..."
            kelly_b, kelly_m, kelly_sig = reglib.run_linmix(x_obs, y_obs, x_err, y_err)

        if method == 'mantz':
            print "Using Mantz Algorithm..."
            mantz_b, mantz_m, mantz_sig = reglib.run_lrgs(x_obs,y_obs,x_err,y_err, nmc=nmc)

    return (kelly_b, kelly_m, kelly_sig), (mantz_b, mantz_m, mantz_sig)

'''
OLD:
def return_methods_list(method):
    #Used to convert a method input choice into a list of used method types.
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
'''

def set_methods_list(method):
    ''' Used to convert a method input choice into a list of used method types.
        Saves to a global variable.
    '''

    global METHODS

    if method is None:
        # Use default method in `param.config`
        method = PARAMETERS['default_methods']

    if method.lower() == 'both':
        METHODS = ['kelly', 'mantz']
    elif method.lower() == 'kelly' or method.lower() == 'mantz':
        # Needed for loop structure
        METHODS = [method]
    else:
        #FIX: If incorrect input, currently uses Kelly as default rather than
        # the one defined in parameter file
        print "WARNNG: Only `kelly`, `mantz`, or `both` are valid method options. Will use Kelly method instead."
        METHODS = ['kelly']

    return


def set_parameters(file):
    ''' Set useful parameters from config file.'''

    global PARAMETERS

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
        name, value = line[0], line[1]

        try:
            # If numeric
            PARAMETERS[name] = eval(value)
        except NameError:
            # If string
            PARAMETERS[name] = value

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
    print('\nChecking dependencies...')
    check_dependencies()

    # Parse all inputted options, regardless of param.config file
    options = parse_opts()

    # Set useful parameters from configure file
    config_file = 'param.config'
    set_parameters(config_file)

    # Set methods list
    method = options.method
    set_methods_list(method)

    # Determine flags
    flags = options.flags

    print('\nInputted options: {}'.format(options))
    print('\nGrabbing data...')

    # Grab and process data from catalog, including flag removal
    data_obs = get_data(options)

    # Scale for linear fitting
    scaled_data = scale(*data_obs)

    print('\nFitting data...')

    # Fit data using linmix, lrgs, or both
    kelly_scaled_fit, mantz_scaled_fit = fit(*scaled_data[:4])
    (x_min, x_max) = (np.min(scaled_data[0]), np.max(scaled_data[0]))

    print('\nMaking plots...')

    # Make all desired plots
    plotlib.make_plots(options, PARAMETERS, METHODS, data_obs, kelly_scaled_fit,
                        mantz_scaled_fit, scaled_data[4], x_min, x_max)

    print('\nDone!')

if __name__ == '__main__':
    main()
