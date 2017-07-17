'''
CluStR takes a FITS cluster catalog and calculates various scaling relations
using either the Kelly method (using `linmix`) or the Mantz method (using
`lrgs`).

Modifies `devon_scaling_relations` by Devon Hollowood.

Spencer Everett, UCSC, 3/2017

Initially based upon `devon_scaling_relations` by Devon Hollowood.
'''

# pylint: disable=invalid-name
# pylint: disable=no-member

import argparse
import os
import cPickle as pickle
import astropy.io.fits as fits
import numpy as np
import reglib  # Regression library
import plotlib  # Plotting library

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
        'lambda': 'lambda',
        'l500kpc': '500_kiloparsecs_band_lumin',
        'lr2500': 'r2500_band_lumin',
        'lr500': 'r500_band_lumin',
        'lr500cc': 'r500_core_cropped_band_lumin',
        't500kpc': '500_kiloparsecs_temperature',
        'tr2500': 'r2500_temperature',
        'tr500': 'r500_temperature',
        'tr500cc': 'r500_core_cropped_temperature'
    }
    return labels[axis_name]


def check_flag(flag):
    ''' Checks whether the flag type is boolean, cutoff, or a range. '''
    # FIX: Flags are now more versatile, so change up the cutoff and

    boolean = {
        'analyzed',
        'detected',
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
        'redshift',
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

    if flag.lower() in boolean:
        try:
            bool_type = PARAMETERS[flag+'_bool_type']
        except KeyError:
            raise TypeError('`{}` is an invalid type entry. Ignoring flag.'
                            .format(flag))
        if bool_type is True or bool_type is False:
            return 'bool'
        else:
            raise TypeError('`{}` is an invalid type entry. Ignoring flag.'
                            .format(flag))

    elif flag.lower() in cut_or_range:
        # Tries cutoff, goes to range if fails. Can only pick one in config!
        try:
            cut_type = PARAMETERS[flag+'_cut_type'].lower()
        except KeyError:
            try:
                range_type = PARAMETERS[flag+'_range_type']
            except KeyError:
                raise TypeError('`{}` is an invalid type entry. Ignoring flag.'
                                .format(flag))
            if range_type == 'inside' or range_type == 'outside':
                return 'range'
            else:
                raise TypeError('`{}` is an invalid type entry. Ignoring flag.'
                                .format(flag))
        if cut_type == 'above' or cut_type == 'below':
            return 'cutoff'
        else:
            raise TypeError('`{}` is an invalid type entry. Ignoring flag.'
                            .format(flag))

    raise TypeError('`{}` is an invalid flag entry. Ignoring flag.'
                    .format(flag))


def create_cuts(data, flags):  # pylint: disable=too-many-branches
    '''
    Apply cuts from `flags` to `data`. Return mask with `True` for all data
    which should be cut
    '''
    mask = np.zeros(len(data), dtype=bool)
    for flag in flags:
        try:
            flag_type = check_flag(flag)
        except TypeError as err:
            print 'WARNING: {}'.format(err)
            print (
                'WARNING: Inputted flag `{}` is not valid. '
                'Ignoring flag for analysis.'
                .format(flag)
            )
            continue

        if flag_type == 'bool':
            bool_type = PARAMETERS[flag + '_bool_type']
            if isinstance(bool_type, bool):
                # keep elements with the boolean value given in `flag`
                cut = data[flag] == (not bool_type)
            else:
                print (
                    'WARNING: Boolean type must be `true` or `false` - '
                    'you entered `{}`. Ignoring `{}` flag.'
                    .format(bool_type, flag)
                )
                continue

        elif flag_type == 'cutoff':
            # Remove data above/below desired cutoff
            cutoff = PARAMETERS[flag + '_cut']
            cut_type = PARAMETERS[flag + '_cut_type'].lower()
            if cut_type == 'above':
                cut = data[flag] < cutoff
            elif cut_type == 'below':
                cut = data[flag] > cutoff
            else:
                print (
                    'WARNING: Cutoff type must be `above` or `below` - '
                    'you entered `{}`. Ignoring `{}` flag.'
                    .format(cut_type, flag)
                )
                continue

        elif flag_type == 'range':
            # Remove data outside of inputted range
            fmin = PARAMETERS[flag + '_range_min']
            fmax = PARAMETERS[flag + '_range_max']
            range_type = PARAMETERS[flag + '_range_type']
            if range_type == 'inside':
                # cut data outside of range
                cut = (data[flag] < fmin) | (data[flag] > fmax)
            elif range_type == 'outside':
                # cut data inside of range
                cut = (data[flag] > fmin) & (data[flag] < fmax)
            else:
                print (
                    'WARNING: Range type must be `inside` or `outside` - '
                    'you entered `{}`. Ignoring `{}` flag.'
                    .format(range_type, flag)
                )
                continue
        mask |= cut
        print (
            'Removed {} clusters due to `{}` flag of type `{}`'
            .format(np.size(np.where(cut)), flag, flag_type)
        )

    return mask


def get_data(options):
    ''' Get x, y, x errors, y errors '''

    hdulist = fits.open(options.catalog)
    data = hdulist[1].data

    label_x = fits_label(options.x)
    label_y = fits_label(options.y)
    x = data[label_x]
    y = data[label_y]

    # Number of original data
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
        # FIX: Should write method to ensure all the counts are what we expect

        mask = create_cuts(data, flags)
        x[mask] = -1
        y[mask] = -1

        print (
            'NOTE: `Removed` counts may be redundant, '
            'as some data fail multiple flags.'
        )

    # Take rows with good data, and all flagged data removed
    good_rows = np.all([x != -1, y != -1], axis=0)

    x = x[good_rows]
    y = y[good_rows]
    x_err = x_err[good_rows]
    y_err = y_err[good_rows]

    print 'Accepted {} data out of {}'.format(np.size(x), N)

    if np.size(x) == 0:
        print (
            '\nWARNING: No data survived flag removal. '
            'Suggest changing flag parameters in `param.config`.'
            '\n\nClosing program...\n'
        )
        raise SystemExit(2)

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


def fit(x_obs, y_obs, x_err, y_err, nmc=5000):
    '''
    Calculates fit using the Kelly (linmix) and/or Mantz (lrgs) methods and
    returns their respective markov chains for the intercept, slope, and sigma.
    Does only Kelly method by default.
    '''

    print "median log(x-x_piv) =", np.median(x_obs)
    # assert np.median(x_obs) == 0.0

    # Set default parameter markov chains to None when both methods aren't used
    kelly_b, kelly_m, kelly_sig = None, None, None
    mantz_b, mantz_m, mantz_sig = None, None, None

    # Iterate through desired methods
    for method in METHODS:
        if method == 'kelly':
            print "Using Kelly Algorithm..."
            kelly_b, kelly_m, kelly_sig = reglib.run_linmix(
                x_obs, y_obs, x_err, y_err
            )

        if method == 'mantz':
            print "Using Mantz Algorithm..."
            mantz_b, mantz_m, mantz_sig = reglib.run_lrgs(
                x_obs, y_obs, x_err, y_err, nmc=nmc
            )

    return (kelly_b, kelly_m, kelly_sig), (mantz_b, mantz_m, mantz_sig)


def set_methods_list(method):
    ''' Used to convert a method input choice into a list of used method types.
        Saves to a global variable.
    '''

    # pylint: disable=global-statement
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
        # FIX: If incorrect input, currently uses Kelly as default rather than
        # the one defined in parameter file
        print (
            'WARNNG: Only `kelly`, `mantz`, or `both` are valid method '
            'options. Will use Kelly method instead.'
        )
        METHODS = ['kelly']

    return


def set_parameters(filename):
    ''' Set useful parameters from config file.'''

    # pylint: disable=global-statement
    global PARAMETERS

    with open(filename) as config_file:
        for line in config_file:
            # Ignore empty lines and comments:
            if line[0:2] == '\n':
                continue
            if line[0] == '#':
                continue
            line.strip()
            line = line.split('#')[0]
            # Remove whitespace and interpret Name:Value pairs:
            line = ''.join(line.split())
            line = line.split(':')
            name, value = line[0], unicode(line[1])

            if value.lower() == 'true':
                PARAMETERS[name] = True
            elif value.lower() == 'false':
                PARAMETERS[name] = False
            elif value.isdigit():
                PARAMETERS[name] = int(value)
            else:
                try:
                    PARAMETERS[name] = float(value)
                except ValueError:
                    # Is string
                    PARAMETERS[name] = value

    return


def check_prefix(options):
    '''
    If no prefix is given, use default set in `param.config`. If default is
    None, then no prefix is given.
    '''
    print 'prefix = {}'.format(options.prefix)
    if options.prefix is None:
        if PARAMETERS['default_prefix'] is None:
            options.prefix = ''
        else:
            options.prefix = PARAMETERS['default_prefix']
    print 'prefix = {}'.format(options.prefix)

    return


def parse_opts():
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser()
    # Required argument for catalog
    parser.add_argument('catalog', help='FITS catalog to open')
    # Required arguement for axes
    valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500',
                  'tr500', 'tr500cc', 'lambda']
    parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
    parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
    # Optional argument for file prefix
    parser.add_argument('-p', '--prefix', help='prefix for output file')
    # Optional argument for regression method(s)
    methods = ['kelly', 'mantz', 'both']
    parser.add_argument(
        '-m',
        '--method',
        help='Choose the `kelly` or `mantz` regression method (or `both`)',
        choices=methods
    )
    # Optional arguments for any flag cuts
    # FIX: in the future, make an allowed choices vector work!
    parser.add_argument(
        '-f',
        '--flags',
        nargs='+',
        type=str,
        help=(
            'Input any desired flag cuts as a list of flag names '
            '(with "" and no spaces!)'
        )
    )
    # Optional argument for which files to be saved
    # FIX: Implement!

    return parser.parse_args()


def check_dependencies():
    '''
    In the future, this function will check if all required packages are
    installed, and, if not, ask if these packages should be downloaded.
    Otherwise exists program.
    '''
    # FIX: Implement!
    # linmix
    # lrgs
    # corner
    # pypdf2
    # check for others!
    pass


def save_data(options, parameters, methods, data_obs, kelly_scaled_fit,
              mantz_scaled_fit, piv, x_min, x_max):
    '''
    Save data locally to a pickle file. Uses default naming scheme if not
    specified in param.config.
    '''
    # pylint: disable=too-many-arguments
    # FIX: refactor this to not have so many arguments

    try: # do try-except instead of if-exists to avoid race condition
        os.makedirs('pickles')
    except os.error: # already existed
        pass

    if parameters['output_filename'] is not None:
        filename = 'pickles/{}'.format(parameters['output_filename'])
        # make sure there is the correct extension
        if filename[-2:] != '.p':
            filename = filename + '.p'
    else:
        filename = 'pickles/Data-{}{}-{}.p'.format(
            options.prefix, fits_label(options.y), fits_label(options.x)
        )

    pickle.dump(
        [options, parameters, methods, data_obs, kelly_scaled_fit,
         mantz_scaled_fit, piv, x_min, x_max],
        open(filename, 'wb')
    )

    return


def main():  # pylint: disable=missing-docstring
    print '\nChecking dependencies...'
    check_dependencies()

    # Parse all inputted options, regardless of param.config file
    options = parse_opts()

    # Set useful parameters from configure file
    config_file = 'param.config'
    set_parameters(config_file)

    # Set default prefix if none entered
    check_prefix(options)

    # Set methods list
    method = options.method
    set_methods_list(method)

    print '\nInputted options: {}'.format(options)
    print '\nGrabbing data...'

    # Grab and process data from catalog, including flag removal
    data_obs = get_data(options)

    # Scale for linear fitting
    scaled_data = scale(*data_obs)

    print '\nFitting data...'

    # Fit data using linmix, lrgs, or both
    kelly_scaled_fit, mantz_scaled_fit = fit(*scaled_data[:4])
    (x_min, x_max) = (np.min(scaled_data[0]), np.max(scaled_data[0]))

    print '\nMaking plots...'

    # Make all desired plots
    plotlib.make_plots(
        options, PARAMETERS, METHODS, data_obs, kelly_scaled_fit,
        mantz_scaled_fit, scaled_data[4], x_min, x_max
    )

    if PARAMETERS['save_data'] is True:
        print '\nSaving data...'
        save_data(options, PARAMETERS, METHODS, data_obs, kelly_scaled_fit,
                  mantz_scaled_fit, scaled_data[4], x_min, x_max)

    print '\nDone!'


if __name__ == '__main__':
    main()
