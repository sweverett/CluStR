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
parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
parser.add_argument('config_file',
    help = 'the filename of the config to run')

#parser.add_argument('plotting_filename',
#    type = str,
#    help = 'the filename of the plotting file to run')
# Optional argument for file prefix
parser.add_argument('-p', '--prefix', help='prefix for output file')
# Optional arguments for any flag cuts
# FIX: in the future, make an allowed choices vector work!
parser.add_argument(
    '-f',
    '--flags',
    nargs='+',
    type=str,
    default=None,
    help=('Input any desired flag cuts as a list of flag names '
    '(with "" and no spaces!)')
)

def fits_label(axis_name):
    ''' Get the FITS column label for `axis_name` '''
#I believe tesla said she just wants one label with no short names
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
        self.flags = args.flags
        self.prefix = args.prefix

        with open(self.filename, 'r') as stream:
            self._config = yaml.safe_load(stream)

        return

    # The following are so we can access the config
    # values similarly to a dict
    def __getitem__(self, key):
        return self._config.get(key)

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
    #take data frome the Catalog Class and pick rows and columns we want to fit
    #dont call flag in main call it here
    def __init__(self, config, catalog):
        self._load_data(config, catalog)

        return

    def _load_data(self, config, catalog):
        '''
        Obtains x, y, x errors, and y errors from config & catalog files.
        '''
        self.xlabel = fits_label(config.x)
        self.ylabel = fits_label(config.y)
        x = catalog[self.xlabel]
        y = catalog[self.ylabel]

        # Size of original data
        N = np.size(x)
        assert N == np.size(y)

        # Scale data if a luminosity
        if config['scale_x_by_ez']:
            x /= Ez(catalog['Redshift'])
        else:
            if self.xlabel[0] == 'l' and self.xlabel != 'lambda':
                print('WARNING: looks like you may be passing a luminosity without'+
                        'setting `scale_x_by_ez: True`. Is that correct?')
        if config['scale_y_by_ez']:
            y /= Ez(catalog['Redshift'])
        else:
            if self.ylabel[0] == 'l' and self.ylabel != 'lambda':
                print('WARNING: looks like you may be passing a luminosity without'+
                        'setting `scale_y_by_ez: True`. Is that correct?')

        self.x_err = (catalog[self.xlabel+'_err_low'] + catalog[self.xlabel+'_err_high']) / 2.
        self.y_err = (catalog[self.ylabel+'_err_low'] + catalog[self.ylabel+'_err_high']) / 2.

        flags = config.flags
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
        x_err = self.x_err[good_rows]
        y_err = self.y_err[good_rows]

        # Cut out any NaNs
        cuts = np.where( (~np.isnan(x)) &
                         (~np.isnan(y)) &
                         (~np.isnan(x_err)) &
                         (~np.isnan(y_err)) )

        self.x = x[cuts]
        self.y = y[cuts]
        self.x_err = x_err[cuts]
        self.y_err = y_err[cuts]

        print('Accepted {} data out of {}'.format(np.size(x), N))

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

        return

    # Plotting the x and y data.
    def plot_data(self, show=True):
        '''
        Plots x-axis and y-axis data.
        '''
        plt.errorbar(self.x, self.y, self.y_err)
        plt.title('{} vs. {}'.format(self.xlabel, self.ylabel))
        plt.xlabel('{}'.format(self.xlabel))
        plt.ylabel('{}'.format(self.ylabel))
        plt.xscale('log')
        plt.yscale('log')
        plt.grid() #add (True, which='both') to display major and minor grids

        if show is True:
            plt.show()

        return

class Fitter:
    def __init__(self):
        self.algorithm = 'linmix'

        return

    def fit(self, data):
        '''
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma.
        '''

        log_x, log_y, log_x_err, log_y_err, piv = self.scale_data(data)

        # run linmix
        print("Using Kelly Algorithm...")
        kelly_b, kelly_m, kelly_sig = reglib.run_linmix(log_x, log_y, log_x_err, log_y_err)

        return kelly_b, kelly_m, kelly_sig

    def scale_data(self, data, piv_type='median'):
        ''' Scale data for fitting'''

        # Log-x before pivot
        xlog = np.log(data.x)

        # Set pivot
        if piv_type == 'median':
            piv = np.log(np.median(data.x))

        log_x = xlog - piv
        log_y = np.log(data.y)

        log_x_err = data.x_err / data.x
        log_y_err = data.y_err / data.y

        return log_x, log_y, log_x_err, log_y_err, piv

"""class SaveData(Fitter):
    def __init__(self, run_options, parameters)
"""

def plot_data(x, y, x_err=None, y_err=None,
              xlabel=None, ylabel=None, log=False, show=True):
    '''
    Independent plotter removed from the Data class
    '''
    plt.plot(x, y, 'o', alpha=0.8, color='tab:blue')
    plt.title('{} vs. {}'.format(xlabel, ylabel))
    plt.xlabel('{}'.format(xlabel))
    plt.ylabel('{}'.format(ylabel))
    if log is True:
        plt.xscale('log')
        plt.yscale('log')
    plt.grid() #add (True, which='both') to display major and minor grids

    if show is True:
        plt.show()

    return

def main():

    args = parser.parse_args()

    config = Config(args)

    catalog = Catalog(args.cat_filename, config)

    data = Data(config, catalog)

    fitter = Fitter()

    b, m, sigma = fitter.fit(data)

    # Scatter plot
    logx, logy, log_x_err, log_y_err, piv = fitter.scale_data(data)
    #logx, logy = fitter.scale_data(data)[0:2]
    #plot_data(logx, logy, show=False, xlabel=config.x, ylabel=config.y)

    xmin = np.min(logx)
    xmax = np.max(logx)
    #dx = abs(xmin - xmax) / 100
    #xx = np.arange(xmin, xmax + dx, dx)
    #plt.plot(xx, np.mean(m)*xx + np.mean(b), lw=3, ls='--', c='k')
    #plt.show()

    print('Done!')

    print('\nFitting Data...')

    plotlib.make_plots(args, config, b, m, sigma, xmax, xmin, piv) #pass on Data class x, y, xerr, yerr?

    return

if __name__ == '__main__':
    main()
