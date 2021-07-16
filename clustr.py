from argparse import ArgumentParser
import os
from astropy.table import Table
import numpy as np
import reglib  # Regression library
import matplotlib.pyplot as plt
import linmix
import yaml
import plotlib
import pyfiglet as pfig
from datetime import datetime


''' Parse command line arguments '''
parser = ArgumentParser()
# Required argument for catalog
parser.add_argument('cat_filename', help='FITS catalog to open')
# Required arguement for axes
valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500',
              'tr500', 'tr500cc', 'lambda', 'lx', 'lam', 'txmm', 'tr2500matcha', 'tr500matcha']
parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
parser.add_argument('config_file',
    help = 'the filename of the config to run')
# Optional argument for file prefix
parser.add_argument('-p', '--prefix', help='prefix for output file')

#----------------------CluStR----------------------------------------

def Ez(z):
    Om = 0.3
    H_0 = 0.7
    h = H_0/100
    return np.sqrt(Om*(1.+z)**3 + h)

# We'll define useful classes here
class Config:
    '''
    Used for CluStR config processing
    Some options:
    - scale_luminosity: Divide luminosity columns by E(z)^-3/2
    '''
    _required_keys = []
    _default_run_name = 'clustr'
    def __init__(self, args):
        # We'll save args as class variables
        self.filename = args.config_file
        self.args = args
        self.x = args.x
        self.y = args.y
        self.prefix = args.prefix

        with open(self.filename, 'r') as stream:
            self._config = yaml.safe_load(stream)

        return

    # The following are so we can access the config
    # values similarly to a dict
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
    #read/load the fits table that contains the data
    def __init__(self,cat_file_name,config):
        self.file_name = cat_file_name

        self._load_catalog()

        return

    def _load_catalog(self):
        self._catalog = Table.read(self.file_name)

        # could do other things...

        return

    # The following are so we can access the catalog
    # values similarly to a dict
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
    '''
    This class takes a catalog table, and grabs only the relevant columns
    for the desired fit using the config dictionary.

    config is expected to act like a dictionary
    '''

    def __init__(self, config, catalog):
        self._load_data(config, catalog)

        return

    def create_cuts(self, config, catalog):
            """
            Apply cuts to data.
            """

            mask = np.zeros(len(catalog), dtype=bool)

            # Boolean Flags
            TF = list(config['Bool_Flag'].keys())[0]

            values = config['Bool_Flag']
            ToF = list(values)[0]
            if ToF != False:
                for bflag_ in config['Bool_Flag']:
                        bool_type = config['Bool_Flag'][bflag_]
                        booleans = list(bool_type.keys())
                        N = len(booleans)
                        a = 0
                        BTF = list(bool_type.values())
                        while(a < N):
                            if isinstance(BTF[a], bool):
                                bflag = booleans[a].replace("_bool_type", "")
                                cutb = catalog[bflag] == (BTF[a])
                                a = a + 1

                            else:
                                print(
                                "Warning: Boolean type must be `True` or  `False` - "
                                "you entered `{}`. Ignoring `{}` flag."
                                .format(bool_type, bflag)
                                )
                                mask |= cutb
                            print(
                            'Removed {} clusters due to `{}` flag of `{}`'
                            .format(np.size(np.where(cutb)), bflag, type(BTF[a-1]))
                            )

            # Cutoff Flags
            for cflag_ in config['Cutoff_Flag']:

                TFc = config['Cutoff_Flag'][cflag_]

                if cflag_ not in ('Other') and list(TFc.keys())[0] != False:
                    cvalues = list(TFc[True].values())
                    cutoff = cvalues[0]
                    cut_type = cvalues[1]

                    if cut_type == 'above':

                        # Nan's interfere with evaluation
                        nan_cut = np.where(np.isnan(catalog[cflag_]))
                        catalog[cflag_][nan_cut] = 2*(cutoff)

                        cutc = catalog[cflag_] < cutoff

                    elif cut_type == 'below':

                        nan_cut = np.where(np.isnan(catalog[cflag_]))
                        catalog[cflag_][nan_cut] = -2*(cutoff)

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
            for rflag_ in config['Range_Flag']:
                TF = config['Range_Flag'][rflag_]
                if rflag_ not in ('Other') and list(TF.keys())[0] != False:

                    rflag = TF[True]

                    for _, rvalues in rflag.items():
                        minmax_ = list(rvalues.values())

                        rmin = minmax_[0]
                        rmax = minmax_[1]
                        range_type = minmax_[2]

                        if range_type == 'inside':
                            cutr = (catalog[rflag_] < rmin) | (catalog[rflag_] > rmax)

                        elif range_type == 'outside':
                            cutr = (catalog[rflag_] > rmin) & (catalog[rflag_] < rmax)

                        else:
                            print (
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
        '''
        Obtains x, y, x errors, and y errors from config & catalog files.
        '''

        x_arg = config.x
        y_arg = config.y
        self.xlabel = config['Column_Names'][x_arg]
        self.ylabel = config['Column_Names'][y_arg]
        x = catalog[self.xlabel]
        y = catalog[self.ylabel]

        # Size of original data
        N = np.size(x)
        assert N == np.size(y)

        # Scale data if a luminosity
        if config['scale_x_by_ez'] == True:
            redshift = config['Redshift']
            x /= Ez(catalog[redshift])
        if config['scale_y_by_ez'] == True:
            redshift = config['Redshift']
            y /= Ez(catalog[redshift])

        # Error Labels
        xlabel_error_low = config["xlabel_err_low"]
        xlabel_error_high = config["xlabel_err_high"]
        ylabel_error_low = config["ylabel_err_low"]
        ylabel_error_high = config["ylabel_err_high"]
        x_err_low = catalog[xlabel_error_low]
        x_err_high = catalog[xlabel_error_high]
        y_err_low = catalog[ylabel_error_low]
        y_err_high = catalog[ylabel_error_high]


        x_err = (catalog[xlabel_error_low] + catalog[xlabel_error_high]) / 2.
        y_err = (catalog[ylabel_error_high] + catalog[ylabel_error_low]) / 2.

        mask = self.create_cuts(config, catalog)

        x[mask] = -1       # All bools_type observations will equal '-1'.

        y[mask] = -1

        print (
        '\nNOTE: `Removed` counts may be redundant, '
        'as some data fail multiple flags.'
        )

        #Take rows with good data, and all flagged data removed
        good_rows = np.all([x != -1, y != -1], axis=0)

        x = x[good_rows]
        y = y[good_rows]
        x_err = x_err[good_rows]
        y_err = y_err[good_rows]
        x_err_low = x_err_low[good_rows]
        x_err_high = x_err_high[good_rows]
        y_err_low = y_err_low[good_rows]
        y_err_high = y_err_high[good_rows]

        # Cut out any NaNs
        cuts = np.where( (~np.isnan(x)) &
                         (~np.isnan(y)) &
                         (~np.isnan(x_err)) &
                         (~np.isnan(y_err)) &
                         (~np.isnan(x_err_low)) &
                         (~np.isnan(x_err_high)) &
                         (~np.isnan(y_err_low)) &
                         (~np.isnan(y_err_high))
                         )
        print(
            'Removed {} nans'
            .format(N - len(cuts[0]))
        )

        self.x = x[cuts]
        self.y = y[cuts]
        self.x_err = x_err[cuts]
        self.y_err = y_err[cuts]
        self.x_err_low = x_err_low[cuts]
        self.x_err_high = x_err_high[cuts]
        self.y_err_low = y_err_low[cuts]
        self.y_err_high = y_err_high[cuts]

        print('Accepted {} data out of {}\n'.format(np.size(self.x), N))

        if np.size(self.x) == 0:
            print (
                '\nWARNING: No data survived flag removal. '
                'Suggest changing flag parameters in `param.config`.'
                '\n\nClosing program...\n'
            )
            raise SystemExit(2)

        #if config.vb is True:
        if config["asymmetric_err"]:
            print('Mean {} error low:'.format(self.xlabel), np.mean(self.x_err_low))
            print('Mean {} error high:'.format(self.xlabel), np.mean(self.x_err_high))
            print('Mean {} error low:'.format(self.ylabel), np.mean(self.y_err_low))
            print('Mean {} error high:'.format(self.ylabel), np.mean(self.y_err_high))

        else:
            print('Mean {} error:'.format(self.xlabel), np.mean(self.x_err))
            print('Mean {} error:'.format(self.ylabel), np.mean(self.y_err))
            print ('\n')

        return

class Fitter:
    """Runs linmix"""

    def __init__(self, data, config):
        self.algorithm = 'linmix'
        self.data_x = data.x
        self.data_y = data.y
        self.data_x_err_obs = data.x_err
        self.data_y_err_obs = data.y_err
        self.data_x_err_low_obs = data.x_err_low
        self.data_x_err_high_obs = data.x_err_high
        self.data_y_err_low_obs = data.y_err_low
        self.data_y_err_high_obs = data.y_err_high
        self.data_xlabel = data.xlabel
        self.data_ylabel = data.ylabel
        self.log_data(data, config)
        self.fit(data)
        self.scaled_fit_to_data(data)

        return

    def fit(self, data):
        '''
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma.
        '''

        # run linmix
        self.kelly_b, self.kelly_m, self.kelly_sigsqr = reglib.run_linmix(
                                                        self.log_x,
                                                        self.log_y,
                                                        self.log_x_err,
                                                        self.log_y_err)

        return

    def log_data(self, data, config):
        ''' Scale data to log'''

        # Log-x before pivot
        xlog = np.log(data.x)

        # Set pivot
        piv_type = config["piv_type"]
        if piv_type == "median":
            self.piv = np.log(np.median(data.x))
        else:
            self.piv = np.log(config['piv_value'])

        else: 
            self.piv = np.log(config["piv_value"])

        self.log_x = xlog - self.piv
        self.log_y = np.log(data.y)

        self.xmin = np.min(self.log_x)
        self.xmax = np.max(self.log_x)

        self.log_x_err = np.log(data.x_err + data.x) - np.log(data.x)
        self.log_y_err = np.log(data.y_err + data.y) - np.log(data.y)

        return

    def scaled_fit_to_data(self, data):
        ''' Get a data set from a scaled fit '''

        #Scale for line fitting
        scaled_x = np.linspace(.5*self.xmin, 1.5*self.xmax, len(self.log_x))
        scaled_y = np.mean(self.kelly_b) + np.mean(self.kelly_m) * scaled_x
        scaled_x_errs = np.zeros(len(self.log_x))
        scaled_y_errs = np.ones(len(self.log_y))*np.mean(self.kelly_m)

        self.unscaled_data = self.unscale(scaled_x, scaled_y, scaled_x_errs, scaled_y_errs, self.piv)

        return

    def unscale(self, x, y, x_err, y_err, x_piv):
        ''' Recover original data from fit-scaled data '''
        ux = np.exp(x + x_piv)
        uy = np.exp(y)
        ux_err = x_err * x
        uy_err = y_err * y

        return (ux, uy, ux_err, uy_err)

    def _recoverXY(self, xObs, xPiv, yObs):
        x = np.exp(xObs+xPiv)
        # x = xObs+xPiv
        y = np.exp(yObs)
        return x, y

    def _regressionLine(self, x, intercept, slope, low, high):
        y = []
        _x = np.linspace(.5*self.xmin, 1.5*self.xmax, len(self.log_x))
        for i, s in zip(self.kelly_b, self.kelly_m):
            y += [i + s * _x]

        y = np.array(y)
        yMed = np.percentile(y, 50, axis=0)
        yLow = np.percentile(y, low, axis=0)
        yUp  = np.percentile(y, high, axis=0)

        return yMed, yUp, yLow

    def _regressionLine_with_scatter(self, low, high):
        y = []
        _x = np.linspace(.5*self.xmin, 1.5*self.xmax, len(self.log_x))
        for i, s, sig in zip(self.kelly_b, self.kelly_m, np.sqrt(self.kelly_sigsqr)):
            y += [i + s * _x + np.random.normal(0.0, sig)]

        y = np.array(y)
        yMed = np.percentile(y, 50, axis=0)
        yLow = np.percentile(y, low, axis=0)
        yUp = np.percentile(y, high, axis=0)
        # print (yMed-yLow)[::5]
        # print (yUp-yMed)[::5]
        return yMed, yUp, yLow

def main():

    #CluStR Banner
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

    #CluStR args
    args = parser.parse_args()

    config = Config(args)

    catalog = Catalog(args.cat_filename, config)

    data = Data(config, catalog)

    fitter = Fitter(data, config)

    print("x-pivot = {}".format(fitter.piv))
    print (
        'Mean b, m, sigsqr: {}, {}, {}'
        .format(np.mean(fitter.kelly_b), np.mean(fitter.kelly_m), np.mean(fitter.kelly_sigsqr))
    )
    print('\n')

    print("Using Kelly Algorithm...")

    print('\nMaking Plots...')

    plotlib.make_plots(args, config, fitter)

    print('Done!')

    return

if __name__ == '__main__':
    main()
