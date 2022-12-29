from argparse import ArgumentParser
from astropy.table import Table
import numpy as np
import reglib  # Regression library
import luminlib
import yaml
import plotlib_2lines
import pyfiglet as pfig
from datetime import datetime
from numpy import savetxt

''' Parse command line arguments '''
parser = ArgumentParser()
# Required argument for catalog
parser.add_argument('cat_filename', help='FITS catalog to open')
# Required argument for axes
valid_axes = ['l500kpc','lr2500s','lx52', 'TR2500', 'LAMBDA', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500', 'tr2500scaled', 'tr500scaled',
              'tr500', 'tr500cc', 'lambda', 'lambdaxmm', 'lambdamatcha', 'lx', 'LAMBDA',
              'lam', 'txmm', 'tr2500matcha', 'tr500matcha', 'tr2500xmm', 'tr500xmm', 'kt', 'lambdachisq', 'R2500',
              'sigma_bi', 'lumin_no_tx',
              'tr500', 'tr500cc', 'lambda', 'lambdaxray', 'lambdachisqxray', 'lambdaxmm', 'lambdamatcha', 'lx',
              'LAMBDA', 'lam', 'txmm', 'tr2500matcha', 'tr500matcha', 'tr2500xmm', 'tr500xmm', 'kt', 'lambdachisq',
              'R2500', 'txmm', 'tmatcha', 'lum', 'txr500matcha', 'txr500xmm']
parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
parser.add_argument('config_file', help='the filename of the config to run')
# Optional argument for file prefix
parser.add_argument('-p', '--prefix', help='prefix for output file')

# ----------------------CluStR----------------------------------------

def Ez(z):
    Om = 0.3
    Ov = 0.7
    return np.sqrt(Om*(1.+z)**3 + Ov)

# We'll define useful classes here
class Config:
    """
    Used for CluStR config processing
    """

    def __init__(self, args):
        """Opens configuration file."""
        self.filename = args.config_file
        self.args = args
        self.x = args.x
        self.y = args.y
        self.prefix = args.prefix

        with open(self.filename, 'r') as stream:
            self._config = yaml.safe_load(stream)

        return

    # Methods used to access values/keys from config.
    def __getitem__(self, key):
        return self._config[key]

    def __setitem__(self, key, value):
        self._config[key] = value

    def __delitem__(self, key):
        del self._config[key]

    def __contains__(self, key):
        return key in self._config

    def __len__(self):
        return len(self._config)

    def __repr__(self):
        return repr(self._config)

class Catalog:
    """
    Read/Load the fits' table that contains the data.
    """

    def __init__(self, cat_file_name):
        self.file_name = cat_file_name

        self._load_catalog()

        return

    def _load_catalog(self):

        """Method used to open catalog."""

        self._catalog = Table.read(self.file_name)

        return

    # Methods used to access values/keys.
    def __getitem__(self, key):
        return self._catalog[key]

    def __setitem__(self, key, value):
        self._catalog[key] = value

    def __delitem__(self, key):
        del self._catalog[key]

    def __contains__(self, key):
        return key in self._catalog

    def __len__(self):
        return len(self._catalog)

    def __repr__(self):
        return repr(self._catalog)

class Data:
    """
    This class takes a catalog table and grabs only the relevant columns
    for the desired fit using the config dictionary.

    Config is expected to act like a dictionary
    """

    def __init__(self, config, catalog):
        self._load_data(config, catalog)

        return

    def create_cuts(self, config, catalog):
        """
        Apply cuts to data. Will remove flags of type Boolean, Cutoff, and Range.
        """

        # Initialize an array of the same size as catalog. Elements are boolean type.
        mask = np.zeros(len(catalog), dtype=bool)

        # Boolean Flags

        # Access True or False key value.
        tf = list(config['Bool_Flag'].keys())[0]

        # Check if user wants boolean cuts.
        if tf:
            # Loop over all boolean flags.
            for bflag_ in config['Bool_Flag'][True]:
                bool_type = config['Bool_Flag'][True][bflag_]

                # Double check if flag is boolean.
                if isinstance(bool_type, bool):

                    bflag = bflag_.replace("_bool_type", "")
                    cutb = catalog[bflag] == bool_type
                else:
                    bflag = bflag_.replace("_bool_type", "")
                    print(
                        "Warning: Boolean type must be `True` or  `False` - "
                        "you entered `{}`. Ignoring `{}` flag."
                        .format(bool_type, bflag)
                    )

                # Include flag cut into mask array.
                mask |= cutb

                print(
                    'Removed {} clusters due to `{}` flag of type boolean.'
                    .format(np.size(np.where(cutb)), bflag_)
                )

        # Cutoff Flags

        # Loop through all cutoffs.
        for cflag_ in config['Cutoff_Flag']:

            TFc = config['Cutoff_Flag'][cflag_]

            # Check if user wants cuts.
            if cflag_ not in 'Other' and list(TFc.keys())[0] is not False:

                # Save values in a list.
                cvalues = list(TFc[True].values())

                # Save in individual variables'
                cutoff = cvalues[0]
                cut_type = cvalues[1]

                # Remove rows below cutoff value.
                if cut_type == 'above':

                    # Nan's interfere with evaluation. Set them to dummy value.
                    nan_cut = np.where(np.isnan(catalog[cflag_]))
                    catalog[cflag_][nan_cut] = 2 * cutoff

                    cutc = catalog[cflag_] < cutoff

                # Remove rows above cutoff value.
                elif cut_type == 'below':

                    # NaN's interfere with evaluation. Set them to dummy value.
                    nan_cut = np.where(np.isnan(catalog[cflag_]))
                    catalog[cflag_][nan_cut] = -2*cutoff

                    cutc = catalog[cflag_] > cutoff

                else:
                    print(
                        'WARNING: Cutoff type must be `above` or `below` - '
                        'you entered `{}`. Ignoring `{}` flag.'
                        .format(cut_type, cflag_))

                mask |= cutc

                print(
                    'Removed {} clusters due to `{}` flag of type cutoff.'
                    .format(np.size(np.where(cutc)), cflag_)
                )

        # Range Flags

        # Loop through all ranges.
        for rflag_ in config['Range_Flag']:
            TF = config['Range_Flag'][rflag_]

            # Check if user wants range cuts.
            if rflag_ not in 'Other' and list(TF.keys())[0] is not False:

                rflag = TF[True]

                for _, rvalues in rflag.items():

                    # Save values to list.
                    minmax_ = list(rvalues.values())

                    rmin = minmax_[0]
                    rmax = minmax_[1]
                    range_type = minmax_[2]

                    # Remove rows outside range.
                    if range_type == 'inside':
                        cutr = (catalog[rflag_] < rmin) | (catalog[rflag_] > rmax)

                    # Remove rows inside range.
                    elif range_type == 'outside':
                        cutr = (catalog[rflag_] > rmin) & (catalog[rflag_] < rmax)

                    else:
                        print(
                            'WARNING: Range type must be `inside` or `outside` - '
                            'you entered `{}`. Ignoring `{}` flag.'
                            .format(range_type, rflag)
                        )
                        continue

                    mask |= cutr

                    print(
                        'Removed {} clusters due to `{}` flag of type range.'
                        .format(np.size(np.where(cutr)), rflag_)
                    )

        return mask

    def _load_data(self, config, catalog):
        """
        Obtains x, y, x errors, and y errors from config & catalog files.
        """

        x_arg = config.x
        y_arg = config.y
        self.xlabel = config['Column_Names'][x_arg]
        self.ylabel = config['Column_Names'][y_arg]


        # Error Labels
        xlabel_error_low = config["xlabel_err_low"]
        xlabel_error_high = config["xlabel_err_high"]
        ylabel_error_low = config["ylabel_err_low"]
        ylabel_error_high = config["ylabel_err_high"]
        x_err_low = catalog[xlabel_error_low]
        x_err_high = catalog[xlabel_error_high]
        y_err_low = catalog[ylabel_error_low]
        y_err_high = catalog[ylabel_error_high]

        x = ((catalog[self.xlabel]))
        y = ((catalog[self.ylabel]))

        # Size of original data
        N = np.size(x)
        assert N == np.size(y)

        # Censored Data
        cenTF = list(config["Censored"].keys())[0]

        if cenTF:
            cenName = config["Censored"][True]
            delta_ = catalog[cenName].astype(np.int64)
            #case to identify type of data
            case = catalog["case"].astype(np.int64)
        else:
            delta_ = np.ones(N)

            case = np.ones(N)
            #case = catalog["case"].astype(np.int64)

        # Cut out any NaNs
        cuts = np.where((~np.isnan(x)) &
                        (~np.isnan(y)) &
                        (~np.isnan(x_err_low)) &
                        (~np.isnan(x_err_high)) &
                        (~np.isnan(y_err_low)) &
                        (~np.isnan(y_err_high))
                        )
        print(
            'Removed {} NaNs'
            .format(N - (N-len(x[cuts])))
        )

        x = x[cuts]
        y = y[cuts]
        x_err_low = x_err_low[cuts]
        x_err_high = x_err_high[cuts]
        y_err_low = y_err_low[cuts]
        y_err_high = y_err_high[cuts]
        delta_ = delta_[cuts]
        case = case[cuts]

        # Scale data
        if config['scale_x_by_ez']:
            redshift = config['Redshift']
            x *= Ez(catalog[redshift][cuts])**(config['scaling_factor_x'][0] / config['scaling_factor_x'][1])

        if config['scale_y_by_ez']:
            redshift = config['Redshift']
            y *= (Ez(catalog[redshift][cuts])**(config['scaling_factor_y']))

        # Set all masked values to negative one.
        mask = self.create_cuts(config, catalog)
        mask = mask[cuts]

        x[mask] = -1

        y[mask] = -1

        print(
            '\nNOTE: `Removed` counts may be redundant, '
            'as some data fail multiple flags.'
        )

        #method for implementing censored data not currently used
        if config['detectedWithTemp']:
            x, y = luminlib.detectedWithTemp(catalog, x=x, y=y)
        if config['detectedWithNoTemp']:
            x, y = luminlib.detectedWithNoTemp(catalog, x=x, y=y, xerr=x_err, yerr=y_err)
        if config['uppLimUndetected']:
            x, y = luminlib.uppLimUndetected(catalog, x=x, y=y)
        else:
            pass

        # Keep rows with good data and remove all flagged data
        good_rows = np.all([x != -1, y != -1], axis=0)

        self.x = x[good_rows]
        self.y = y[good_rows]
        self.x_err_low = x_err_low[good_rows]
        self.x_err_high = x_err_high[good_rows]
        self.y_err_low = y_err_low[good_rows]
        self.y_err_high = y_err_high[good_rows]
        self.delta_ = delta_[good_rows]
        self.case = case[good_rows]

        #saving columns for later
        #np.savetxt('casesredlow.csv', self.case, delimiter=',')


        print('Accepted {} data out of {}\n'.format(np.size(self.x), N))

        if np.size(self.x) == 0:
            print(
                '\nWARNING: No data survived flag removal. '
                'Suggest changing flag parameters in `param.config`.'
                '\n\nClosing program...\n'
            )
            raise SystemExit(2)

        # if config is True:
        print(f'Mean {self.xlabel} error low: {np.mean(self.x_err_low)}')
        print(f'Mean {self.xlabel} error high: {np.mean(self.x_err_high)}')
        print(f'Mean {self.ylabel} error low: {np.mean(self.y_err_low)}')
        print(f'Mean {self.ylabel} error high: {np.mean(self.y_err_high)}')


        return

class Fitter:
    """Runs linmix algorithm using the regression library."""

    def __init__(self, data, config):
        """ Here we can use the super method to inherit
            the attributes from the Data class.
        """

        self.scaled_x = None
        self.log_y_err = None
        self.log_x_err = None
        self.xmax = None
        self.xmin = None
        self.log_y = None
        self.log_x = None
        self.piv = None
        self.mean_sigsqr = None
        self.mean_slope = None
        self.mean_int = None
        self.kelly_sigsqr = None
        self.kelly_m = None
        self.kelly_b = None
        self.algorithm = 'linmix'
        self.data_x = data.x
        self.data_y = data.y
        self.data_x_err_low_obs = data.x_err_low
        self.data_x_err_high_obs = data.x_err_high
        self.data_y_err_low_obs = data.y_err_low
        self.data_y_err_high_obs = data.y_err_high
        self.data_xlabel = data.xlabel
        self.data_ylabel = data.ylabel
        self._constant = config['scale_line']
        self.case_n = data.case
        self.log_data(config)
        self.fit(data)
        self.scaled_fit_to_data()
        return

    def fit(self, data):
        """
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma_sqr.
        """

        self.kelly_b, self.kelly_m, self.kelly_sigsqr = reglib.run_linmix(
                                                            x=self.log_x,
                                                            y=self.log_y,
                                                            err_x=self.log_x_err,
                                                            err_y=self.log_y_err,
                                                            delta=data.delta_)

#saving to add to plotlib when plotting 2 lines
        #savetxt('kelly_bT_r2500_joint_0.2_to_0.4.csv', self.kelly_b, delimiter=',')
        #savetxt('kelly_mT_r2500_joint_0.2_to_0.4.csv', self.kelly_m, delimiter=',')
        #savetxt('kelly_sigsqrT_r2500_joint_0.2_to_0.4.csv', self.kelly_sigsqr, delimiter=',')

        self.mean_int = np.mean(self.kelly_b)
        self.mean_slope = np.mean(self.kelly_m)
        self.mean_sigsqr = np.mean(self.kelly_sigsqr)
        return

    def log_data(self, config):
        """ Scale data to log"""
        # Set pivot
        piv_type = config["piv_type"]
        if piv_type == "median":
            self.piv = np.log(np.median(self.data_x))
        else:
            self.piv = np.log(config['piv_value'])

        #find symmetric errors

        log_y_max = np.log(self.data_y + self.data_y_err_high_obs)
        log_y_min = np.log(self.data_y - self.data_y_err_low_obs)

        log_x_max = np.log(self.data_x + self.data_x_err_high_obs) - self.piv
        log_x_min = np.log(self.data_x - self.data_x_err_low_obs) - self.piv

        #centralize x and y and divide x by log pivot

        self.log_y = (log_y_max + log_y_min)/2
        self.log_x = (log_x_max + log_x_min)/2

        self.log_y_err = log_y_max - self.log_y
        self.log_x_err = log_x_max - self.log_x

        self.xmin = np.min(self.log_x)
        self.xmax = np.max(self.log_x)

        self.xlim = [0.07*np.min(self.data_x), 1.3*np.max(self.data_x)]
        self.xPlot = np.linspace(np.log(self.xlim[0])-self.piv, np.log(self.xlim[1])-self.piv, 201)

        return

    def scaled_fit_to_data(self):
        """ Calculate scaled linear values. """
        #self.scaled_x = np.linspace(np.log(self.xlim[0])-self.piv, np.log(self.xlim[1])-self.piv, 98)

        #if conf interval is not extending all the way adjust these
        self.scaled_x = np.linspace(1.7*self.xmin, 1.5*self.xmax, len(self.log_x))
        scaled_y = self.mean_int + self.mean_slope * self.scaled_x
        scaled_x_errs = np.zeros(len(self.log_x))
        scaled_y_errs = np.ones(len(self.log_y))*self.mean_slope

        return self.scaled_x, scaled_y, scaled_x_errs, scaled_y_errs

    def unscaled(self):
        """ Recover original data from scaled_fit_to_data() """

        # Grab log-scaled linear values.
        sx, sy, sx_err, sy_err = self.scaled_fit_to_data()

        # Recover
        ux = np.exp(sx + self.piv)
        uy = np.exp(sy)
        ux_err = sx_err * sx
        uy_err = sy_err * sy

        return ux, uy, ux_err, uy_err

    def _recoverY(self, yObs):
        """This method will return unscaled Y."""
        y = np.exp(yObs)
        return y

    def confInterval(self, low, high):
        """This method will calculate confidence interval from y distribution."""

        y = []
        _x = np.linspace(.02*self.xmin, 10000*self.xmax, len(self.log_x))
        for i, s in zip(self.kelly_b, self.kelly_m):
            y += [i + s * self.scaled_x]

        y = np.array(y)
        yMed = np.percentile(y, 50, axis=0)
        yLow = np.percentile(y, low, axis=0)
        yUp = np.percentile(y, high, axis=0)

        return yMed, yLow, yUp

    def sigmaBands(self, low, high):
        """ This method calculates sigma bands."""

        y = []
        _x = np.linspace(.02*self.xmin, 10000*self.xmax, len(self.log_x))
        for i, s, sig in zip(self.kelly_b, self.kelly_m, (self.kelly_sigsqr)):
            y += [i + s * self.scaled_x + np.random.normal(0.0, sig)]

        y = np.array(y)
        yMed = np.percentile(y, 50, axis=0)
        yLow = np.percentile(y, low, axis=0)
        yUp = np.percentile(y, high, axis=0)

        return yMed, yLow, yUp

class Banner:
    """Contains Program Banner"""

    def __init__(self):
        # CluStR Banner
        ascii_banner = pfig.figlet_format("CluStR")
        print(ascii_banner)
        print("-----------------------------------")
        print("This package calculates various \nscaling relations from cluster catalogs.")
        print("\n")

        # Returns the current local date
        now = datetime.now()
        print(now)
        print("-----------------------------------")
        print("\n")

def main():

    # CluStR Banner
    Banner()

    # CluStR args
    args = parser.parse_args()

    config = Config(args)

    catalog = Catalog(args.cat_filename)

    data = Data(config, catalog)

    fitter = Fitter(data, config)
    print(f"x-pivot = {fitter.piv}")
    print(f"Mean Intercept: {np.mean(fitter.kelly_b)}")
    print(f"Mean Slope: {np.mean(fitter.kelly_m)}")
    print(f"Mean Variance: {np.mean(fitter.kelly_sigsqr)}")

    print('\n')

    print("Using Kelly Algorithm...")

    print('\nMaking Plots...')

    plotlib_2lines.make_plots(args, config, fitter)

    print('Done!')

    return


if __name__ == '__main__':
    main()
