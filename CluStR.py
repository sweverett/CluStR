'''
CluStR takes a FITS cluster catalog and calculates various scaling relations
using either the Kelly method (using `linmix`) or the Mantz method (using `lrgs`).
Modifies `devon_scaling_relations` by Devon Hollowood.

Spencer Everett, UCSC, 3/2017
'''

#pylint: disable=invalid-name
#pylint: disable=no-member

import argparse
import matplotlib, os, corner, ast
matplotlib.use('Agg')
import matplotlib.pylab as plt
import astropy.io.fits as fits
import numpy as np
#import eslib
#import lrgslib
import reglib
import PyPDF2

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
        'offset_500',
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
        print('WARNING: {} is an invalid flag entry. Ignoring entry.'.format(flag))

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
                pass
            elif flag_type == 'range':
                # Remove data outside of inputted range
                pass
            else:
                print('WARNING: Inputted flag {} is not valid. Ignoring flag for analysis.'.format(flag))
                continue

            '''
            #OLD
            flagged = np.where(data[flag] == True)
            x[flagged] = -1
            y[flagged] = -1
            L = np.size(flagged)
            print('Removed {} data points due to {} flag'.format(L,flag))
            '''

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

def axis_label(axis_name):
    ''' Get plot axis name for `axis_name` '''
    labels = {
        'lambda':'richness',
        'l500kpc':'500 Kiloparsec Soft-Band Luminosity (ergs/s)',
        'lr2500':'r2500 Soft-Band Luminosity (ergs/s)',
        'lr500':'r500 Soft-Band Luminosity (ergs/s)',
        'lr500cc':'Core-Cropped r500 Soft-Band Luminosity (ergs/s)',
        't500kpc':'500 Kiloparsec Temperature (keV)',
        'tr2500':'r2500 Temperature (keV)',
        'tr500':'r500 Temperature (keV)',
        'tr500cc':'Core-Cropped r500 Temperature (keV)'
    }
    return labels[axis_name]

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def plot_scatter(options, data_obs, kelly_scaled_fit, mantz_scaled_fit,x_piv, x_min, x_max):
    ''' Plot data '''
    (x_obs, y_obs, x_err_obs, y_err_obs) = data_obs

    # Plot data
    plt.errorbar(x_obs, y_obs, xerr=x_err_obs, yerr=y_err_obs, c='r', fmt='o')

    methods = return_methods_list(options.method)

    for method in methods:
        if method == 'kelly':
            scaled_fit = kelly_scaled_fit
            color = 'b'

        elif method == 'mantz':
            scaled_fit = mantz_scaled_fit
            color = 'm'

        (fit_int, fit_slope, fit_sig) = scaled_fit
        data_fit = scaled_fit_to_data(x_min, x_max, x_piv, scaled_fit)
        (x_fit, y_fit, _, _) = data_fit

        # plot fit
        plt.loglog(
            x_fit, y_fit, color, linewidth=2.0,
            label=(r'$({0:0.2g} \pm {1:0.2g}) (x/x_{{piv}})^{{{2:0.2f} \pm {3:0.2f}}} '
                r'(\sigma = {4:0.2f} \pm {5:0.2f})$'
            ).format(
                np.exp(np.mean(fit_int)),
                np.exp(np.mean(fit_int)) * np.std(fit_int),
                np.mean(fit_slope),
                np.std(fit_slope),
                np.mean(fit_sig),
                np.std(fit_sig)
            )
        )
        '''
        FIX: Re-implement, and default to having no name attached for method
        plt.loglog(
            x_fit, y_fit, color, linewidth=2.0,
            label=('{0}:'
                r'$({1:0.2g} \pm {2:0.2g}) (x/x_{{piv}})^{{{3:0.2g} \pm {4:0.2g}}} '
                r'(\sigma^2 = {5:0.2g} \pm {6:0.2g})$'
            ).format(
                method.capitalize(),
                np.exp(np.mean(fit_int)),
                np.exp(np.mean(fit_int)) * np.std(fit_int),
                np.mean(fit_slope),
                np.std(fit_slope),
                np.mean(fit_sig),
                np.std(fit_sig)
            )
        )
        '''

    plt.xlabel(axis_label(options.x), fontsize=16)
    plt.ylabel(axis_label(options.y), fontsize=16)
    plt.xlim([0.8*np.min(x_obs), 1.2*np.max(x_obs)])
    plt.ylim([0.8*np.min(y_obs), 1.2*np.max(y_obs)])

    plt.legend(loc=2)

    plt.savefig('Scatter-{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)),
                bbox_inches='tight')

    return

def plot_corners(options, kelly_scaled_fit, mantz_scaled_fit,burn=0):
    '''
    Makes corner plots for the desired Kelly and/or Mantz method parameter posteriors. Burn is the
    burn in period parameter.
    '''

    # FIX: add a condition to adapt burn # if needed

    methods = return_methods_list(options.method)

    N = np.size(methods) # Number of subplots
    n = 1 # Subplot counter

    #figures = []

    for method in methods:
        if method == 'kelly':
            scaled_fit = kelly_scaled_fit
        elif method == 'mantz':
            scaled_fit = mantz_scaled_fit
        else:
            print "WARNNG: Only `kelly`, `mantz`, or `both` are valid method options. Will use Kelly method instead."
            scaled_fit = kelly_scaled_fit

        # Set up subplot
        plt.subplot(N, 1, n)

        (B, M, S) = scaled_fit

        # Paramter Limits
        blo, bhi = min(B[burn:]), max(B[burn:])
        mlo, mhi = min(M[burn:]), max(M[burn:])
        slo, shi = min(S[burn:]), max(S[burn:])

        #FIX: maybe use lo = -hi for symmetry?? Can cause issues for small min

        sf = 0.25 # scale factor
        db = sf*abs(bhi - blo)
        dm = sf*abs(mhi - mlo)
        ds = sf*abs(shi - slo)
        blo, bhi = blo-db, bhi+db
        mlo, mhi = mlo-db, mhi+dm
        slo, shi = slo-ds, shi+ds
        #blo, bhi = blo-abs(blo)*sf, bhi+abs(bhi)*sf
        #mlo, mhi = mlo-abs(mlo)*sf, mhi+abs(mhi)*sf
        #slo, shi = slo-abs(slo)*sf, shi+abs(shi)*sf


        data = np.transpose( (B, M, S) )

        fig = corner.corner(data, labels=['b','m','s'], range=[(blo,bhi),(mlo,mhi),(slo,shi)],quantiles=[0.16,0.5,0.84],
                        show_titles=True, title_args={"fontsize": 18},
                        plot_datapoints=True, fill_contours=True, levels=[0.68, 0.95], color='b', bins=40, smooth=1.0);

        fig.suptitle('{} Method'.format(method.capitalize()),fontsize=16)

        plt.savefig('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))#,bbox_inches='tight')

        #figures.append(fig)

        n += 1 # Iterate counter

    #plt.gcf().set_size_inches(20,10)

    return

def plot_chains(options,kelly_scaled_fit,mantz_scaled_fit,burn=0):
    '''
    Use this to examine chain convergence. May implement convergence tests in future.
    '''

    methods = return_methods_list(options.method)

    for method in methods:
        if method == 'kelly':
            scaled_fit = kelly_scaled_fit
        elif method == 'mantz':
            scaled_fit = mantz_scaled_fit
        else:
            print "WARNNG: Only `kelly`, `mantz`, or `both` are valid method options. Will use Kelly method instead."
            scaled_fit = kelly_scaled_fit

        # Initialize
        B, M, S = None, None, None

        # Unpack fit parameters
        (B, M, S) = scaled_fit
        # Remove burn-in period
        B, M, S = B[burn:], M[burn:], S[burn:]
        # Take averages
        b, m, s = np.mean(B), np.mean(M), np.mean(S)

        # Length of chain
        nmc = np.size(B)

        fig = plt.figure()

        plt.subplot(311)

        plt.plot(M,'o',markerfacecolor="None")
        #plt.plot((0,nmc), (1.0, 1.0), 'k--', linewidth=2)
        plt.plot((0,nmc),(m,m),'r--')
        plt.xlabel('Chain Number')
        plt.ylabel('Slope')
        #plt.gcf().set_size_inches(8,4)

        plt.subplot(312)

        plt.plot(B,'o',markerfacecolor="None")
        #plt.plot((0,nmc), (0.0, 0.0), 'k--', linewidth=2)
        plt.plot((0,nmc),(b,b),'r--')
        plt.xlabel('Chain Number')
        plt.ylabel('Intercept')
        #plt.gcf().set_size_inches(8,4)

        plt.subplot(313)

        plt.plot(S,'o',markerfacecolor="None")
        #plt.plot((0,nmc), (9.0, 9.0), 'k--', linewidth=2)
        plt.plot((0,nmc),(s,s),'r--')
        plt.xlabel('Chain Number')
        plt.ylabel(r'$\sigma^2$')
        #plt.gcf().set_size_inches(8,4)

        fig.suptitle('Markov Chains for {} Method'.format(method.capitalize()),fontsize=16)

        fig.set_size_inches(10,10)
        #plt.show()

        plt.savefig('Chains-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

        plt.clf()

    return

def make_plots(options, data_obs, kelly_scaled_fit, mantz_scaled_fit, piv, x_min, x_max):
    '''Calls both plotting functions and then combines all outputs into a single PDF.'''

    # FIX: Re-implement, and allow plotting choices as arguments!

    plot_scatter(options, data_obs, kelly_scaled_fit, mantz_scaled_fit,piv, x_min, x_max)
    plot_corners(options,kelly_scaled_fit,mantz_scaled_fit)
    plot_chains(options,kelly_scaled_fit,mantz_scaled_fit)

    # Initialize pdf Luminosity
    pdfs = []

    # Add scatter/fit plot
    pdfs.append('Scatter-{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    # Add corner plots
    methods = return_methods_list(options.method)

    for method in methods:
        pdfs.append('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    # Separate loops for the right order
    for method in methods:
        pdfs.append('Chains-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    merger = PyPDF2.PdfFileMerger()

    for pdf in pdfs:
        merger.append(pdf)

    merger.write('{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    # Unless otherwise specified, delete individual plots
    ## FIX: for now, always delete
    if True:
        for pdf in pdfs:
            os.remove(pdf)

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
    # linmix
    # lrgs
    # corner
    # pypdf2
    pass

def main(): #pylint: disable=missing-docstring
    print('Checking dependencies...')
    check_dependencies()
    options = parse_opts()
    method = options.method
    flags = options.flags
    print('Inputted options:',options)
    print('Grabbing data...')
    data_obs = get_data(options)
    scaled_data = scale(*data_obs)
    print('Fitting data...')
    kelly_scaled_fit, mantz_scaled_fit = fit(method,*scaled_data[:4])
    (x_min, x_max) = (np.min(scaled_data[0]), np.max(scaled_data[0]))
    print('Making plots...')
    make_plots(options, data_obs, kelly_scaled_fit, mantz_scaled_fit, scaled_data[4], x_min, x_max)
    #plot_chains('kelly',kelly_scaled_fit)
    #plot_chains('mantz',mantz_scaled_fit)
    print('Done!')

if __name__ == '__main__':
    main()
