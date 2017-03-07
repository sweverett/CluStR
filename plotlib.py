# Plotting library for CluStR

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt


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

    # OLD: now replaced by METHODS
    #methods = clustr.return_methods_list(options.method)

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

        # plot fit
        if parameters['show_method_name'] is True or len(METHODS) > 1:
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

def plot_corners(options, kelly_scaled_fit, mantz_scaled_fit,burn=0):
    '''
    Makes corner plots for the desired Kelly and/or Mantz method parameter posteriors. Burn is the
    burn in period parameter.
    '''

    # FIX: add a condition to adapt burn # if needed

    # OLD: now uses METHODS
    #methods = clustr.return_methods_list(options.method)

    N = np.size(methods) # Number of subplots
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

        fig.suptitle('{} Method'.format(method.capitalize()),fontsize=16)

        plt.savefig('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))#,bbox_inches='tight')

        n += 1 # Iterate counter

    return

def plot_chains(options,kelly_scaled_fit,mantz_scaled_fit,burn=0):
    '''
    Use this to examine chain convergence. May implement convergence tests in future.
    '''

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

def make_plots(options, data_obs, kelly_scaled_fit, mantz_scaled_fit, piv, x_min, x_max):
    '''Calls both plotting functions and then combines all outputs into a single PDF.'''

    # OLD: now uses METHODS
    # Retreive methods list
    #methods = clustr.return_methods_list(options.method)

    if parameters['scatter'] is True:
        plot_scatter(options, data_obs, kelly_scaled_fit, mantz_scaled_fit,piv, x_min, x_max)
        # Add scatter/fit plot
        pdfs.append('Scatter-{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    if parameters['corner'] is True:
        plot_corners(options,kelly_scaled_fit,mantz_scaled_fit)
        # Add corner plot(s)
        for method in METHODS:
            pdfs.append('Corner-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    if parameters['chains'] is True:
        plot_chains(options,kelly_scaled_fit,mantz_scaled_fit)
        # Add chain plot(s)
        for method in METHODS:
            pdfs.append('Chains-{}-{}{}-{}.pdf'.format(method,options.prefix, fits_label(options.y), fits_label(options.x)))

    if parameters['residuals'] is True:
        # FIX: Implement!
        pass

    # Initialize pdf list
    pdfs = []
    merger = PyPDF2.PdfFileMerger()

    for pdf in pdfs:
        merger.append(pdf)

    # Save combined output file
    merger.write('{}{}-{}.pdf'.format(options.prefix, fits_label(options.y), fits_label(options.x)))

    # Unless otherwise specified, delete individual plots
    if parameters['save_all_plots'] is True:
        for pdf in pdfs:
            os.remove(pdf)
