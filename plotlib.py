"""Plotting library for CluStR """

import corner
import PyPDF2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.ticker import LogFormatter, ScalarFormatter, FuncFormatter
# import seaborn as ssb
# plt.style.use('seaborn')
# matplotlib.use('Agg')

# pylint: disable=invalid-name

# ----------------------------------------------------------------------
# Individual plotting functions

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

def plot_scatter(args, fitter, config):
    """ Plot data """

    # Grab data references
    # Symmetric Errors
    x_obs = fitter.data_x
    y_obs = fitter.data_y
    x_err_low_obs = fitter.data_x_err_low_obs
    x_err_high_obs = fitter.data_x_err_high_obs
    y_err_low_obs = fitter.data_y_err_low_obs
    y_err_high_obs = fitter.data_y_err_high_obs

    # Asymmetric Errors
    # x_err_obs_low = fitter.data_x_err_low_obs
    # x_err_obs_high = fitter.data_x_err_high_obs
    # y_err_obs_low = fitter.data_y_err_low_obs
    # y_err_obs_high = fitter.data_y_err_high_obs

    # x_err_obs_asym = [x_err_obs_low, x_err_obs_high]
    # y_err_obs_asym = [y_err_obs_low, y_err_obs_high]

    # Plot data
    fig, ax = plt.subplots()
    plt.errorbar(x_obs, y_obs, xerr=np.array([x_err_low_obs, x_err_high_obs]), yerr=np.array([y_err_low_obs, y_err_high_obs]),
                 ecolor='k',
                 fmt='bo',
                 lw=1,
                 markersize=2,
                 markeredgecolor='k',
                 capsize=1,
                 label='_nolegend_'
                 )
    print('Reporting Symmetric Error Bars.')

    # Grab linmix data
    fit_int, fit_slope, fit_sig = fitter.kelly_b, fitter.kelly_m, fitter.kelly_sigsqr

    # Line data
    (x_fit, y_fit, _, _) = fitter.unscaled()

    # Plot Linear Fit (x_fit = unscaled x) and (y_fit = unscaled line)
    plt.loglog(
        x_fit, y_fit, base=np.e, color='navy', linewidth=2.0,
        label=(
            r'$({0:0.2g} \pm {1:0.2g})'
            r'(x/x_{{piv}})^{{{2:0.2f} \pm {3:0.2f}}}'
            r'(\sigma^2 = {4:0.2f} \pm {5:0.2f})$'
        ).format(
            np.exp(np.mean(fit_int)),
            np.exp(np.mean(fit_int)) * np.std(fit_int),
            np.mean(fit_slope),
            np.std(fit_slope),
            np.mean(fit_sig),
            np.std(fit_sig)
        )
    )
            
    # Confidence Interval
    yMed0, yLow0, yUp0 = fitter.confInterval(16, 84)
    yMed0 = fitter._recoverY(yMed0)
    yUp0 = fitter._recoverY(yUp0)
    yLow0 = fitter._recoverY(yLow0)
    plt.fill_between(x_fit, yUp0, yLow0, color='b', alpha=0.3, label=r'68% Confidence Interval')

    # Sigma Bands
    yMed1, yLow1, yUp1 = fitter.sigmaBands(16, 84)
    yUp1 = fitter._recoverY((yUp1 - yMed1) + yMed1)
    yLow1 = fitter._recoverY((yLow1 - yMed1) + yMed1)
    plt.fill_between(x_fit, yUp1, yLow1, color='teal', alpha=0.25, label= r'1$\sigma$ Band')

    yMed2, yLow2, yUp2 = fitter.sigmaBands(16, 84)
    yUp2 = fitter._recoverY(2*(yUp2 - yMed2) + yMed2)
    yLow2 = fitter._recoverY(2*(yLow2 - yMed2) + yMed2)
    plt.fill_between(x_fit, yUp2, yLow2, color='teal', alpha=0.2, label= r'2$\sigma$ Band')

    # -----------------------------------------------------------------
    # Plot Labels
    if list(config["Plot_Labels"].keys())[0]:
        xname = list(config["Plot_Labels"][True].values())[0]
        yname = list(config["Plot_Labels"][True].values())[1]

    else:
        xname = fitter.data_xlabel.capitalize()
        yname = fitter.data_ylabel

    def myLogFormat(y,f):
        # Find the number of decimal places required
        decimalplaces = int(np.maximum(-np.log10(y/10.0+0.01),0)) # =0  numbers >=1
        # Insert that number into a format string
        formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
        # Return the formatted tick label
        return formatstring.format(y)
    #set labels for x and y axes
    ax.set_xlabel(f'{xname}', fontsize=10)
    ax.set_ylabel(f'{yname}', fontsize=10)
    #set limits for x and y axes
    ax.set_xlim([0.7*np.min(x_obs), 1.3*np.max(x_obs)])
    ax.set_ylim([0.4, 22])
    #set scale for axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    #set tick parameters
    plt.tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='major',      # both major and minor ticks are affected
                    bottom=True,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=True     # labels along the bottom edge are off
                    ) 
    plt.tick_params(
                    axis='y',          
                    which='major',      
                    bottom=True,      
                    top=False,         
                    labelbottom=True
                    )
    
    #set ticks for x and y axis
    ax.set_xticks([20, 40, 60, 80, 100, 200])
    ax.set_yticks([0.6, 0.8, 1, 2, 4, 6, 8, 10, 20])
    
    #Honestly not sure what this does but the ticks don't format correctly without it
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
    ax.grid(which='major', color='k', alpha=0.2)
    ax.grid(which='minor', color='k', alpha=0.1)
    ax.legend(loc='best', fontsize='x-small')

    plt.savefig(
        'Scatter-{}{}-{}.pdf'
        .format(
            args.prefix,
            fitter.data_ylabel,
            fitter.data_xlabel
        ),
        bbox_inches='tight'
    )

    return

def plot_residuals(args, fitter):
    """
    FIX: Description
    """

    (lx_obs, ly_obs, _lx_err_obs, _ly_err_obs) = fitter.log_x, fitter.log_y, fitter.log_x_err, fitter.log_y_err
    _x_piv = fitter.piv

    (B, M, _S) = fitter.kelly_b, fitter.kelly_m, fitter.kelly_sigsqr

    b, m = np.mean(B), np.mean(M)

    # Calculate residuals
    x_fit = lx_obs
    y_fit = m*x_fit+b
    # FIX: Find out which normalization to use!
    # residuals = (ly_obs - y_fit) / ly_err_obs
    residuals = (ly_obs - y_fit) / np.std(ly_obs)

    '''
    # Make residual plot fig1 = plt.figure(1) #Plot Data-model frame1 =
    fig1.add_axes((.1,.3,.8,.6)) #xstart, ystart, xend, yend [units are
    fraction of the image frame, from bottom left corner]
    plt.errorbar(lx_obs, ly_obs, xerr=lx_err_obs, yerr=ly_err_obs, c='r',
    fmt='o') plt.plot(x_fit,y_fit,'b') #Best fit model
    frame1.set_xticklabels([]) #Remove x-tic labels for the first frame

    #Residual plot frame2 = fig1.add_axes((.1,.1,.8,.2))
    plt.plot(lx_obs,residuals,'ob')
    plt.plot([np.min(lx_obs),np.max(lx_obs)],[0,0],'k--',linewidth=2)

    plt.title('{} Method'.format(method.capitalize()),fontsize=14)
    '''

    # Bin number
    # FIX: make bins automatically consistent with Michigan group
    nbin = 18

    plt.style.use('seaborn')
    plt.hist(residuals, nbin)
    plt.xlabel(r'$\Delta(\ln X)/\sigma_{\ln X}$', fontsize=11)
    plt.ylabel('Count', fontsize=11)

    plt.title(
        '{} Residuals'
        .format(fitter.data_ylabel),
        fontsize=11
    )

    plt.savefig(
        'Residuals-{}{}-{}.pdf'
        .format(
            args.prefix,
            fitter.data_ylabel,
            fitter.data_ylabel
        )
    )

    return

def plot_corners(args, config, fitter):
    """
    Makes corner plots for the desired Kelly method parameter
    posteriors. Burn is the burn in period parameter.
    """

    burn = config['burn']

    # FIX: Is this still being used?
    N = np.size(9)  # Number of subplots
    n = 1  # Subplot counter

    # Set up subplot
    plt.style.use('seaborn')
    plt.subplot(N, 1, n)

    (B, M, S) = fitter.kelly_b, fitter.kelly_m, fitter.kelly_sigsqr

    # Parameter Limits
    blo, bhi = min(B[burn:]), max(B[burn:])
    mlo, mhi = min(M[burn:]), max(M[burn:])
    slo, shi = min(S[burn:]), max(S[burn:])

    # FIX: maybe use lo = -hi for symmetry?? Can cause issues for small min

    sf = 0.25  # scale factor
    db = sf*abs(bhi - blo)
    dm = sf*abs(mhi - mlo)
    ds = sf*abs(shi - slo)
    blo, bhi = blo-db, bhi+db
    mlo, mhi = mlo-db, mhi+dm
    slo, shi = slo-ds, shi+ds
    # blo, bhi = blo-abs(blo)*sf, bhi+abs(bhi)*sf
    # mlo, mhi = mlo-abs(mlo)*sf, mhi+abs(mhi)*sf
    # slo, shi = slo-abs(slo)*sf, shi+abs(shi)*sf

    data1 = np.transpose((B, M, S))

    fig = corner.corner(
        data1,
        labels=['b', 'm', 's'],
        range=[(blo, bhi), (mlo, mhi), (slo, shi)],
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_args={"fontsize": 18},
        plot_datapoints=True,
        fill_contours=False,
        levels=[0.68, 0.95],
        color='k',
        bins=40,
        smooth=1.0
    )
    fig.suptitle('Posterior Distribution',
                 fontsize=14)

    plt.savefig(
        'Corner-{}{}-{}.pdf'
        .format(
            args.prefix,
            fitter.data_ylabel,
            fitter.data_xlabel
        )
    )

    n += 1  # Iterate counter

    return


def plot_chains(args, config, fitter):
    """
    Use this to examine chain convergence. May implement convergence tests in
    the future.
    """

    burn = config['burn']

    # Initialize
    B, M, S = None, None, None

    # Unpack fit parameters
    (B, M, S) = fitter.kelly_b, fitter.kelly_m, fitter.kelly_sigsqr
    # Remove burn-in period
    B, M, S = B[burn:], M[burn:], S[burn:]
    # Take averages
    b, m, s = np.mean(B), np.mean(M), np.mean(S)

    # Length of chain
    nmc = np.size(B)

    plt.style.use('ggplot')
    fig = plt.figure()

    plt.subplot(311)

    plt.plot(M, 'o', markerfacecolor="None")
    plt.plot((0, nmc), (m, m), 'r--')
    plt.xlabel('Chain Number')
    plt.ylabel('Slope')

    plt.subplot(312)

    plt.plot(B, 'o', markerfacecolor="None")
    plt.plot((0, nmc), (b, b), 'r--')
    plt.xlabel('Chain Number')
    plt.ylabel('Intercept')

    plt.subplot(313)

    plt.plot(S, 'o', markerfacecolor="None")
    plt.plot((0, nmc), (s, s), 'r--')
    plt.xlabel('Chain Number')
    plt.ylabel(r'$\sigma^2$')

    fig.suptitle(
        '{} vs. {} \n\nMarkov Chains for Kelly Method'
        .format(fitter.data_ylabel,
                fitter.data_xlabel
                ),
        fontsize=16
    )

    fig.set_size_inches(10, 10)

    plt.savefig(
        'Chains-{}{}-{}.pdf'
        .format(
            args.prefix,
            fitter.data_ylabel,
            fitter.data_xlabel
        )
    )

    plt.clf()

    return


# ----------------------------------------------------------------------
# Make all plots

def make_plots(args, config, fitter):
    """
    Calls both plotting functions and then combines all outputs into a single
    PDF.
    """

    # pylint: disable=global-statement

    # OLD: now uses METHODS
    # Retrieve methods list

    # Initialize pdf list
    pdfs = []

    if config['scatter'] is True:
        plot_scatter(args, fitter, config)
        # Add scatter/fit plot
        pdfs.append(
            'Scatter-{}{}-{}.pdf'
            .format(
                args.prefix,
                fitter.data_ylabel,
                fitter.data_xlabel
            )
        )

    if config['residuals'] is True:
        plot_residuals(args, fitter, config)
        # Add residual plot(s)
        pdfs.append(
            'Residuals-{}{}-{}.pdf'.format(
                args.prefix,
                fitter.data_ylabel,
                fitter.data_xlabel
            )
        )

    if config['corner'] is True:
        plot_corners(args, config, fitter)
        # Add corner plot(s)
        pdfs.append(
            'Corner-{}{}-{}.pdf'
            .format(
                args.prefix,
                fitter.data_ylabel,
                fitter.data_xlabel
            )
        )

    if config['chains'] is True:
        plot_chains(args, config, fitter)
        # Add chain plot(s)
        pdfs.append(
            'Chains-{}{}-{}.pdf'
            .format(
                args.prefix,
                fitter.data_ylabel,
                fitter.data_xlabel
            )
        )
    if config['save_all_plots'] is True:
        merger = PyPDF2.PdfFileMerger()

        for pdf in pdfs:
            merger.append(pdf)

        # Save combined output file
        merger.write(
            '{}{}-{}.pdf'
            .format(args.prefix, fitter.data_ylabel, fitter.data_xlabel)
        )
    else:
        pass

    return
