# Plotting library for CluStR

from clustr import fits_label, scale
import matplotlib
matplotlib.use('Agg')
import PyPDF2, corner, os
import numpy as np
import matplotlib.pylab as plt

# ----------------------------------------------------------------------
# Some useful global variables

# Parameter list used throughout. See `param.config`
PARAMETERS = {}

# Method List
METHODS = []

# ----------------------------------------------------------------------
# Helper functions

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

def scaled_fit_to_data(x_min, x_max, x_piv, scaled_fit):
    ''' Get a data set from a scaled fit '''
    (fit_int, fit_slope, fit_sig) = scaled_fit
    scaled_x = np.linspace(x_min, x_max, 101)
    scaled_y = np.mean(fit_int) + np.mean(fit_slope) * scaled_x
    scaled_x_errs = np.zeros(101)
    scaled_y_errs = np.ones(101)*np.mean(fit_sig)
    unscaled_data = unscale(scaled_x, scaled_y, scaled_x_errs, scaled_y_errs, x_piv)
    return unscaled_data

def unscale(x, y, x_err, y_err, x_piv):
    ''' Recover original data from fit-scaled data '''
    return (np.exp(x + x_piv), np.exp(y), x_err * x, y_err * y)

def check_convergence():
    '''FIX: In future, will use autocorrelations to check convergence of MCMC chains.'''
    pass

# ----------------------------------------------------------------------
# Inividual plotting functions

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def plot_scatter(options, data_obs, kelly_scaled_fit, mantz_scaled_fit,x_piv, x_min, x_max):
    ''' Plot data '''
    (x_obs, y_obs, x_err_obs, y_err_obs) = data_obs

    # Plot data
    plt.errorbar(x_obs, y_obs, xerr=x_err_obs, yerr=y_err_obs, c='r', fmt='o')

    for method in METHODS:
        if method == 'kelly':
            scaled_fit = kelly_scaled_fit
            color = 'b'

        elif method == 'mantz':
            scaled_fit = mantz_scaled_fit
            color = 'm'

        (fit_int, fit_slope, fit_sig) = scaled_fit
        data_fit = scaled_fit_to_data(x_min, x_max, x_piv, scaled_fit)
        (x_fit, y_fit, _, _) = data_fit

        print('mean b, m, sig: {}, {}, {}'.format(np.mean(fit_int),np.mean(fit_slope),np.mean(fit_sig)))

        # plot fit
        if PARAMETERS['show_method_name'] is True or len(METHODS) > 1:
            # Prints method name in legend (default if more than 1 method)
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

        else:
            # Doesn't print method label
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

    plt.xlabel(axis_label(options.x), fontsize=16)
    plt.ylabel(axis_label(options.y), fontsize=16)
    plt.xlim([0.8*np.min(x_obs), 1.2*np.max(x_obs)])
    plt.ylim([0.8*np.min(y_obs), 1.2*np.max(y_obs)])

    plt.legend(loc=2)

    plt.savefig('Scatter-{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)),
                bbox_inches='tight')

    return

def plot_residuals(options, data_obs, kelly_scaled_fit, mantz_scaled_fit):
    '''
    FIX: Description
    '''

    #(x_obs, y_obs, x_err_obs, y_err_obs) = data_obs
    (lx_obs, ly_obs, lx_err_obs, ly_err_obs, x_piv) = scale(*data_obs)

    for method in METHODS:
        if method == 'kelly':
            scaled_fit = kelly_scaled_fit
        elif method == 'mantz':
            scaled_fit = mantz_scaled_fit
        else:
            print "WARNNG: Only `kelly`, `mantz`, or `both` are valid method options. Will use Kelly method instead."
            scaled_fit = kelly_scaled_fit

        (B, M, S) = scaled_fit

        b, m, sig = np.mean(B), np.mean(M), np.mean(S)

        # Calculate residuals
        x_fit = lx_obs
        y_fit = m*x_fit+b
        # FIX: Find out which normalization to use!
        #residuals = (ly_obs - y_fit) / ly_err_obs
        residuals = (ly_obs - y_fit) / np.std(ly_obs)

        '''
        # Make residual plot
        fig1 = plt.figure(1)
        #Plot Data-model
        frame1 = fig1.add_axes((.1,.3,.8,.6))
        #xstart, ystart, xend, yend [units are fraction of the image frame, from bottom left corner]
        plt.errorbar(lx_obs, ly_obs, xerr=lx_err_obs, yerr=ly_err_obs, c='r', fmt='o')
        plt.plot(x_fit,y_fit,'b') #Best fit model
        frame1.set_xticklabels([]) #Remove x-tic labels for the first frame

        #Residual plot
        frame2 = fig1.add_axes((.1,.1,.8,.2))
        plt.plot(lx_obs,residuals,'ob')
        plt.plot([np.min(lx_obs),np.max(lx_obs)],[0,0],'k--',linewidth=2)

        plt.title('{} Method'.format(method.capitalize()),fontsize=14)
        '''

        fig = plt.figure()

        # Bin number
        # FIX: make bins automatically consistent with Michigan group
        nbin = 18

        plt.hist(residuals,nbin)
        plt.xlabel(r'$\Delta(\ln X)/\sigma_{\ln X}$',fontsize=16)
        plt.ylabel('Count',fontsize=16)

        if PARAMETERS['show_method_name'] is True or len(METHODS) > 1:
            plt.title('{} Residuals for {} Method'.format(fits_label(options.y),method.capitalize()),fontsize=16)
        else:
            plt.title('{} Residuals'.format(fits_label(options.y)),fontsize=16)


        plt.savefig('Residuals-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    return


def plot_corners(options, kelly_scaled_fit, mantz_scaled_fit):
    '''
    Makes corner plots for the desired Kelly and/or Mantz method parameter posteriors. Burn is the
    burn in period parameter.
    '''

    burn = PARAMETERS['burn']

    # FIX: Is this still being used?
    N = np.size(METHODS) # Number of subplots
    n = 1 # Subplot counter

    for method in METHODS:
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

        if PARAMETERS['show_method_name'] is True or len(METHODS) > 1:
            fig.suptitle('Posterior for {} Method'.format(method.capitalize()),fontsize=16)
        else:
            fig.suptitle('Posterior',fontsize=16)

        plt.savefig('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))#,bbox_inches='tight')

        n += 1 # Iterate counter

    return

def plot_chains(options,kelly_scaled_fit,mantz_scaled_fit):
    '''
    Use this to examine chain convergence. May implement convergence tests in future.
    '''

    burn = PARAMETERS['burn']

    # OLD: now uses MEHTODS
    #methods = clustr.return_methods_list(options.method)

    for method in METHODS:
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

# ----------------------------------------------------------------------
# Make all plots

def make_plots(options, parameters, methods, data_obs, kelly_scaled_fit, mantz_scaled_fit, piv, x_min, x_max):
    '''Calls both plotting functions and then combines all outputs into a single PDF.'''

    # OLD: now uses METHODS
    # Retreive methods list
    #methods = clustr.return_methods_list(options.method)

    # Sets module parameters to those set in clustr.py
    global PARAMETERS
    global METHODS
    PARAMETERS = parameters
    METHODS = methods

    # Initialize pdf list
    pdfs = []

    if PARAMETERS['scatter'] is True:
        plot_scatter(options, data_obs, kelly_scaled_fit, mantz_scaled_fit,piv, x_min, x_max)
        # Add scatter/fit plot
        pdfs.append('Scatter-{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    if PARAMETERS['residuals'] is True:
        plot_residuals(options, data_obs, kelly_scaled_fit, mantz_scaled_fit)
        # Add residual plot(s)
        for method in METHODS:
            pdfs.append('Residuals-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    if PARAMETERS['corner'] is True:
        plot_corners(options,kelly_scaled_fit,mantz_scaled_fit)
        # Add corner plot(s)
        for method in METHODS:
            pdfs.append('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    if PARAMETERS['chains'] is True:
        plot_chains(options,kelly_scaled_fit,mantz_scaled_fit)
        # Add chain plot(s)
        for method in METHODS:
            pdfs.append('Chains-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    merger = PyPDF2.PdfFileMerger()

    for pdf in pdfs:
        merger.append(pdf)

    # Save combined output file
    merger.write('{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    # Unless otherwise specified, delete individual plots
    if parameters['save_all_plots'] is False:
        for pdf in pdfs:
            os.remove(pdf)

    return
