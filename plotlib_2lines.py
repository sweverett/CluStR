x_datax_data"""Plotting library for CluStR """
import corner
import PyPDF2
import csv
from astropy.table import Table, join, vstack
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.ticker import LogFormatter, ScalarFormatter, FuncFormatter


def data1():
    file1 = "joint_target.fits"
    catalog1 = Table.read(file1)
    x_data = catalog1["lambda"]
    y_data = catalog1["r2500_temperature_scaled"]
    y_data_err_low = catalog1["r2500_temperature_err_low_scaled"]
    y_data_err_high = catalog1["r2500_temperature_err_high_scaled"]
    x_data_err_low = catalog1["lambda_err_low"]
    x_data_err_high = catalog1["lambda_err_high"]

    #THIS NEEDS TO BE FIXED TOO MATCH CLUSTR ERR HANDLING

    #change these
    return  x_data, y_data, y_data_err_low, y_data_err_high, x_data_err_low, x_data_err_high

def log_data1():
    """ Scale data to log"""

    # Log-x before pivot
    x_data, y_data, y_data_err_low, y_data_err_high, x_data_err_low, x_data_err_high = data1()
    y_max = y_data + y_data_err_high
    y_min = y_data - y_data_err_low

    x_max = x_data + x_data_err_high
    x_min = x_data- x_data_err_low

    log_y_err_high = np.log(y_max) - np.log(y_data)
    log_y_err_low = np.log(y_data) - np.log(y_min)
    log_x_err_high = np.log(x_max) - np.log(x_data)
    log_x_err_low = np.log(x_data) - np.log(x_min)

    #symmetric log errors
    log_y_err = (log_y_err_high + log_y_err_low)/2
    log_x_err = (log_x_err_high + log_x_err_low)/2

    #centralize x and y

    #and divide x by pivot
    log_y_max = np.log(y_data) + log_y_err_high
    log_y_min = np.log(y_data) - log_y_err_low
    log_x_max = np.log(x_data) + log_x_err_high
    log_x_min = np.log(x_data) - log_x_err_low


    log_y = (log_y_min + log_y_max)/2
    log_x = (log_x_min + log_x_max)/2 - np.log(70)

    xmin = np.min(log_x)
    xmax = np.max(log_x)

    xlim = [0.7*np.min(x_data), 1.3*np.max(x_data)]
    xPlot = np.linspace(np.log(self.xlim[0])-self.piv, np.log(self.xlim[1])-self.piv, 201)

    xlog = np.log(x_data1)

    # Set pivot

    piv = np.log(70)

    log_x = xlog - piv
    log_y = np.log(y_data1)

    xmin = np.min(log_x)
    xmax = np.max(log_y)

    log_x_err = np.log(x_err_obs1 + x_data1) - xlog
    log_y_err = np.log(y_err_obs1 + y_data1) - log_y

    return log_x, log_y, log_x_err, log_y_err, xmin, xmax

def scaled_fit_to_data1():
    """ Calculate scaled linear values. """
    log_x1, log_y1, log_x_err1, log_y_err1, xmin1, xmax1 = log_data1()
    x_data2, y_data2, x_err_obs2, y_err_obs2 = data1()
    xlim = [0.7*np.min(x_data2), 1.3*np.max(x_data2)]
    scaled_x = np.linspace(np.log(xlim[0])-np.log(70), np.log(xlim[1])-np.log(70), 162)
    scaled_y = 1.46 + 0.41 * scaled_x
    scaled_x_errs = np.zeros(len(log_x1))
    scaled_y_errs = np.ones(len(log_y1))*0.41

    return scaled_x, scaled_y, scaled_x_errs, scaled_y_errs

def unscaled1():
    """ Recover original data from scaled_fit_to_data() """

    # Grab log-scaled linear values.
    sx, sy, sx_err, sy_err = scaled_fit_to_data1()

    # Recover to cartesian
    ux = np.exp(sx + np.log(70))
    uy = np.exp(sy)
    ux_err = sx_err * sx
    uy_err = sy_err * sy

    return ux, uy, ux_err, uy_err

def _recoverY1(yObs):
    """This method will return unscaled Y."""
    y = np.exp(yObs)
    return y

kelly_b1 = []
with open("kelly_bTa.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        kelly_b1.append(row)
kelly_m1 = []
with open("kelly_mTa.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        kelly_m1.append(row)
kelly_s1 = []
with open("kelly_sigsqrTa.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        kelly_s1.append(row)



print(np.mean(kelly_b1))
def confInterval1(low, high):
    """This method will calculate confidence interval from y distribution."""
    log_x1, log_y1, log_x_err1, log_y_err1, xmin1, xmax1 = log_data1()
    scaled_x, scaled_y, scaled_x_errs, scaled_y_errs = scaled_fit_to_data1()
    y = []

    #NEED TO RUN LINMIX IN HERE
    for i, s in zip(kelly_b1, kelly_m1):
        y += [i + s * scaled_x]

    y = np.array(y)
    yMed = np.percentile(y, 50, axis=0)
    yLow = np.percentile(y, low, axis=0)
    yUp = np.percentile(y, high, axis=0)

    return yMed, yUp, yLow

def sigmaBands1(low, high):
    """ This method calculates sigma bands."""
    log_x1, log_y1, log_x_err1, log_y_err1, xmin1, xmax1 = log_data1()
    y = []

    scaled_x, scaled_y, scaled_x_errs, scaled_y_errs = scaled_fit_to_data1()
    for i, s, sig in zip(kelly_b1, kelly_m1, (kelly_s1)):
        y += [i + s * scaled_x + np.random.normal(0.0, sig)]

    y = np.array(y)
    yMed = np.percentile(y, 50, axis=0)
    yLow = np.percentile(y, low, axis=0)
    yUp = np.percentile(y, high, axis=0)

    return yMed, yUp, yLow


def plot_scatter(args, fitter, config):
    """ Plot data """

    # Grab data references
    # Symmetric Errors
    case = fitter.case_n
    x_obs = fitter.data_x
    y_obs = fitter.data_y
    x_err_obs = fitter.data_x_err_obs
    y_err_obs = fitter.data_y_err_obs

#THIS NEEDS TO BE CHANGED FROM THIS METHOD TO TAKING THE
    x_obs1, y_obs1, x_err_obs1, y_err_obs1 = data1()

    #categories, Chandra data to be red, XMM black. represent upper limits (delta = 0) with a caret (red+black)
    fig, ax = plt.subplots()

#Case 1 + 2:
    plt.errorbar(x_obs1, y_obs1, xerr=x_err_obs1, yerr=y_err_obs1,
                 ecolor='r',
                 color = 'r',
                 fmt='o',
                 lw=.75,
                 elinewidth = .5,
                 markersize=1,
                 markeredgecolor='r',
                 capsize=.5,
                 label='Target'
                 )


    plt.errorbar(x_obs, y_obs, xerr=x_err_obs, yerr=y_err_obs,
                 ecolor='k',
                 color = 'k',
                 fmt='o',
                 lw=.75,
                 elinewidth = .5,
                 markersize=1,
                 markeredgecolor='k',
                 capsize=.5,
                 label='Serendipitous'
                 )


    # Grab linmix data
    fit_int, fit_slope, fit_sig = fitter.kelly_b, fitter.kelly_m, fitter.kelly_sigsqr





    fit_int1, fit_slope1, fit_sig1 = np.mean(kelly_b1), np.mean(kelly_m1), np.mean(kelly_s1)
    print(fit_int1, fit_slope1, fit_sig1)
    # Line data
    (x_fit, y_fit, _, _) = fitter.unscaled()
    (x_fit1, y_fit1, _, _) = unscaled1()

    plt.loglog(
        x_fit1, y_fit1, basex=np.e, basey=np.e, color='maroon', linewidth=2.0,
        label=(
            r'$Target ({0:0.2g} \pm {1:0.2g})'
            r'(x/x_{{piv}})^{{{2:0.2f} \pm {3:0.2f}}}'
            r'(\sigma^2 = {4:0.2f} \pm {5:0.2f})$'
        ).format(
            np.exp(np.mean(fit_int1)),
            np.exp(np.mean(fit_int1)) * np.std(fit_int1),
            np.mean(fit_slope1),
            np.std(fit_slope1),
            np.mean(fit_sig1),
            np.std(fit_sig1)
        )
    )
    #ignore
        # Confidence Interval
    yMed01, yLow01, yUp01 = confInterval1(16, 84)
    yMed01 = _recoverY1(yMed01)
    yUp01 = _recoverY1(yUp01)
    yLow01 = _recoverY1(yLow01)
    plt.fill_between(x_fit1, yUp01, yLow01, color='b', alpha=0.3, label=r'68% Confidence Interval Target')

        # Sigma Bands
    yMed11, yLow11, yUp11 = sigmaBands1(16, 84)
    yUp11 = _recoverY1((yUp11 - yMed11) + yMed11)
    yLow11 = _recoverY1((yLow11 - yMed11) + yMed11)
    plt.fill_between(x_fit1, yUp11, yLow11, color='dimgrey', alpha=0.25, label= r'1$\sigma$ Band Target')

    yMed21, yLow21, yUp21 = sigmaBands1(16, 84)
    yUp21 = _recoverY1(2*(yUp21 - yMed21) + yMed21)
    yLow21 = _recoverY1(2*(yLow21 - yMed21) + yMed21)
    plt.fill_between(x_fit1, yUp21, yLow21, color='dimgrey', alpha=0.2, label= r'2$\sigma$ Band Target')


    # Plot Linear Fit (x_fit = unscaled x) and (y_fit = unscaled line)
    plt.loglog(
        x_fit, y_fit, basex=np.e, basey=np.e, color='navy', linewidth=2.0,
        label=(
            r'$Serendipitous ({0:0.2g} \pm {1:0.2g})'
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
#ignore
    # Confidence Interval
    yMed0, yLow0, yUp0 = fitter.confInterval(16, 84)
    yMed0 = fitter._recoverY(yMed0)
    yUp0 = fitter._recoverY(yUp0)
    yLow0 = fitter._recoverY(yLow0)
    plt.fill_between(x_fit, yUp0, yLow0, color='b', alpha=0.3, label=r'68% Confidence Interval Serendipitous')

    # Sigma Bands
    yMed1, yLow1, yUp1 = fitter.sigmaBands(16, 84)
    yUp1 = fitter._recoverY((yUp1 - yMed1) + yMed1)
    yLow1 = fitter._recoverY((yLow1 - yMed1) + yMed1)
    plt.fill_between(x_fit, yUp1, yLow1, color='teal', alpha=0.25, label= r'1$\sigma$ Band Serendipitous')

    yMed2, yLow2, yUp2 = fitter.sigmaBands(16, 84)
    yUp2 = fitter._recoverY(2*(yUp2 - yMed2) + yMed2)
    yLow2 = fitter._recoverY(2*(yLow2 - yMed2) + yMed2)
    plt.fill_between(x_fit, yUp2, yLow2, color='teal', alpha=0.2, label= r'2$\sigma$ Band Serendipitous')

    # -----------------------------------------------------------------
    # Plot Labels
#ignore
    if list(config["Plot_Labels"].keys())[0]:
        xname = list(config["Plot_Labels"][True].values())[0]
        yname = list(config["Plot_Labels"][True].values())[1]

    else:
        xname = fitter.data_xlabel.capitalize()
        yname = fitter.data_ylabel
#ignore
    def myLogFormat(y,f):
        # Find the number of decimal places required
        decimalplaces = int(np.maximum(-np.log10(y/10.0+0.01),0)) # =0  numbers >=1
        # Insert that number into a format string
        formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
        # Return the formatted tick label
        return formatstring.format(y)

    ax.set_xlabel(f'{xname}', fontsize=10)
    ax.set_ylabel(f'{yname}', fontsize=10)
    #we will work on this together
    ax.set_xlim([0.7*np.min(x_obs), 1.5*np.max(x_obs)])
    ax.set_ylim([0.5, 22])
    ax.set_ylim([.1*np.min(y_obs), 2.5*np.max(y_obs1)]) # (ignore for now)
    #we will work on this together
    ax.set_xscale('log', subsx=[2, 4, 6, 8,10])
    ax.set_yscale('log', subsy=[2, 4, 6, 8, 10,200])
    ax.set_yticks([20,30])
    ax.set_xticks([20,40,300])
    #ax.tick_params(axis='both', which='major', direction='in', length=8, width=1.)
    #ax.tick_params(axis='both', which='minor', direction='in', length=4, width=0.5)
    ax.xaxis.set_major_formatter(LogFormatter())
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    #ax.yaxis.set_major_formatter(LogFormatter())
    #ax.yaxis.set_minor_formatter(FuncFormatter(myLogFormat))

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
    plt.show()
    return
'''
    plt.errorbar(x_obs[0:74], y_obs[0:74], xerr=x_err_obs[0:74], yerr=y_err_obs[0:74],
                 ecolor='r',
                 color = 'r',
                 fmt='o',
                 lw=.75,
                 elinewidth = .5,
                 markersize=1,
                 markeredgecolor='r',
                 capsize=.5,
                 label='Chandra'
                 )
                #5

    plt.errorbar(x_obs[281:441], y_obs[281:441], xerr=x_err_obs[281:441], yerr=y_err_obs[281:441],
                    ecolor='k',
                    color = 'k',
                    fmt='bo',
                     lw=.75,
                     elinewidth = .5,
                     markersize=1,
                     markeredgecolor='k',
                     capsize=.5,
                     label='XMM'
                     )
#Case 3 + 4:
    plt.errorbar(x_obs[74:281], y_obs[74:281], xerr=x_err_obs[74:281], yerr=y_err_obs[74:281],uplims=True,
                    ecolor='r',
                    color = 'r',
                    fmt='^',
                     lw=.75,
                     elinewidth = .5,
                     markersize=2.8,
                     markeredgecolor='r',
                     capsize=.5,
                     label='Chandra Upper Limits'
                     )



#Case 6 + 7 + 8:
    plt.errorbar(x_obs[441:684], y_obs[441:684], xerr=x_err_obs[441:684], yerr=y_err_obs[441:684], uplims=True,
                    ecolor='k',
                    color = 'k',
                    fmt='^',
                     lw=.75,
                     elinewidth = .5,
                     markersize= 2.8,
                     markeredgecolor='k',
                     capsize= .5,
                     label='XMM Upper Limits'
                     )
'''

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
