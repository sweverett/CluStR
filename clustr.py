from argparse import ArgumentParser
import os
import pickle as cPickle
from astropy.table import Table
import numpy as np
import reglib  # Regression library
import matplotlib.pyplot as plt
import linmix
import yaml
import plotlib

#import pudb

from astropy.io import fits

''' Parse command line arguments '''
parser = ArgumentParser()
# Required argument for catalog
parser.add_argument('cat_filename', help='FITS catalog to open')
# Required arguement for axes
valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500',
              'tr500', 'tr500cc', 'lambda']
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

class Data(Catalog):
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

            maskb = np.zeros(len(catalog), dtype=bool)
            maskc = np.zeros(len(catalog), dtype=bool)
            maskr = np.zeros(len(catalog), dtype=bool)

            # Boolean Flags
            for bflag_ in config['Bool_Flag']:
                bool_type = config['Bool_Flag'][bflag_]

                if isinstance(bool_type, bool):

                    bflag = bflag_.replace("_bool_type", "")

                    cutb = catalog[bflag] == (bool_type)

                else:
                    print(
                        "Warning: Boolean type must be `True` or  `False` - "
                        "you entered `{}`. Ignoring `{}` flag."
                        .format(bool_type, bflag)
                    )

                maskb |= cutb
                print(
                    'Removed {} clusters due to `{}` flag of `{}`'
                    .format(np.size(np.where(cutb)), bflag_, type(bool_type))
                )

            # Cutoff Flags
            for cflag_ in config['Cutoff_Flag']:

                TFc = config['Cutoff_Flag'][cflag_]

                if cflag_ not in ('Other') and TFc.keys()[0] != False:
                    cvalues = TFc[True].values()
                    cut_type = cvalues
                    cutoff = cvalues[0]

                    if cut_type == 'above':

                        # Nan's interfere with evaluation

                        cutc = catalog[cflag_] < cutoff

                    elif cut_type == 'below':

                        cutc = catalog[cflag_] > cutoff

                    else:
                        print(
                            'WARNING: Cutoff type must be `above` or `below` - '
                            'you entered `{}`. Ignoring `{}` flag.'
                            .format(cut_type, cflag_))

                    maskc |= cutc

                    print(
                        'Removed {} clusters due to `{}` flag of `{}`'
                        .format(np.size(np.where(cutc)), cflag_, type(cflag_))
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
                        #print(range_type)

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

                        maskr |= cutr

                        print(
                            'Removed {} clusters due to `{}` flag of `{}`'
                            .format(np.size(np.where(cutr)), rflag_, type(range_type))
                        )

            return maskb, maskc, maskr

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
        if config['scale_x_by_ez']:
            x /= Ez(catalog['Z_1'])     #Changed Redshift == Z_1
        else:
            if self.xlabel[0] == 'L' and self.xlabel != 'LAM':
                print('WARNING: looks like you may be passing a luminosity without'+
                        'setting `scale_x_by_ez: True`. Is that correct?')
        if config['scale_y_by_ez']:
            y /= Ez(catalog['Z_1'])
        else:
            if self.ylabel[0] == 'l' and self.ylabel != 'lambda':
                print('WARNING: looks like you may be passing a luminosity without'+
                        'setting `scale_y_by_ez: True`. Is that correct?')

        self.x_err = (catalog[self.xlabel+'m'] + catalog[self.xlabel+'p']) / 2.     # Changed _err_low == m
        self.y_err = (catalog[self.ylabel+'m'] + catalog[self.ylabel+'p']) / 2.     # Changed _err_high == p

        maskb, maskc, maskr = self.create_cuts(config, catalog)

        x[maskb] = -1       # All bools_type observations will equal '-1'.
        x[maskr] = -1       # All range_type observations will equal '-1'
        x[maskc] = -1       # All cut_type observations will equal '-1'.

        y[maskb] = -1
        y[maskr] = -1
        y[maskc] = -1

        print (
        '\nNOTE: `Removed` counts may be redundant, '
        'as some data fail multiple flags.'
        )

        # Take rows with good data, and all flagged data removed
        good_rows = np.all([x != -1, y != -1], axis=0)

        x = x[good_rows]
        y = y[good_rows]
        x_err = self.x_err[good_rows]
        y_err = self.y_err[good_rows]

        # Cut out any NaNs
        cuts = np.where( (~np.isnan(x)) &
                         (~np.isnan(y)) &
                         (~np.isnan(x_err)) &
                         (~np.isnan(y_err)) )
        print(
            'Removed {} nans'
            .format(np.size(np.where(cuts)))
        )

        self.x = x[cuts]
        self.y = y[cuts]
        self.x_err = x_err[cuts]
        self.y_err = y_err[cuts]

        print('Accepted {} data out of {}\n'.format(np.size(self.x), N))

        if np.size(x) == 0:
            print (
                '\nWARNING: No data survived flag removal. '
                'Suggest changing flag parameters in `param.config`.'
                '\n\nClosing program...\n'
            )
            raise SystemExit(2)

        #if config.vb is True:
        print('mean x error:', np.mean(self.x_err))
        print('mean y error:', np.mean(self.y_err))
        print ('\n')

        return

class Fitter:
    """Runs linmix"""

    def __init__(self):
        self.algorithm = 'linmix'

        return

    def fit(self, data):
        '''
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma.
        '''

        log_x, log_y, log_x_err, log_y_err = self.scale_data(data)[0:4]

        # run linmix
        kelly_b, kelly_m, kelly_sig = reglib.run_linmix(log_x, log_y, log_x_err, log_y_err)

        return kelly_b, kelly_m, kelly_sig

    def scale_data(self, data, piv_type='median'):
        ''' Scale data for fitting'''

        # Log-x before pivot
        xlog = np.log(data.x)

        # Set pivot
        if piv_type == 'median':
            self.piv = np.log(np.median(data.x))

        log_x = xlog - self.piv
        log_y = np.log(data.y)

        xmin = np.min(log_x)
        xmax = np.max(log_x)

        log_x_err = data.x_err / data.x
        log_y_err = data.y_err / data.y

        return log_x, log_y, log_x_err, log_y_err, xmin, xmax, self.piv

def main():

    args = parser.parse_args()

    config = Config(args)

    catalog = Catalog(args.cat_filename, config)

    data = Data(config, catalog)

    fitter = Fitter()

    print("x-pivot = {}".format(fitter.scale_data(data)[6]))
    print('\n')

    print("Using Kelly Algorithm...")

    print('\nMaking Plots...')

    plotlib.make_plots(args, config, data, fitter)

    print('Done!')

    return

if __name__ == '__main__':
    main()
