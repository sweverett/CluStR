from argparse import ArgumentParser
import os
from astropy.table import Table
import numpy as np
import reglib  # Regression library
import luminlib
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
valid_axes = ['l500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500', 'tr2500scaled', 'tr500scaled',
              'tr500', 'tr500cc', 'lambda', 'lambdaxmm', 'lambdamatcha', 'lx', 'LAMBDA', 'LAMBDAXRAY',
              'lam', 'txmm', 'tr2500matcha', 'tr500matcha', 'tr2500xmm', 'tr500xmm', 'kt', 'lambdachisq','R2500', 'sigma_bi', 'lumin_no_tx',
              'tr500', 'tr500cc', 'lambda', 'lambdaxray','lambdachisqxray','lambdaxmm', 'lambdamatcha', 'lx', 'LAMBDA',
              'lam', 'txmm', 'tr2500matcha', 'tr500matcha', 'tr2500xmm', 'tr500xmm', 'kt', 'lambdachisq','R2500'
              'txmm', 'tmatcha', 'lum', 'txr500matcha', 'txr500xmm', 'lum_no_tx', 'LX52', 'lupper']
parser.add_argument('x', help='what to plot on x axis', choices=valid_axes)
parser.add_argument('y', help='what to plot on y axis', choices=valid_axes)
parser.add_argument('config_file',
    help = 'the filename of the config to run')
# Optional argument for file prefix
parser.add_argument('-p', '--prefix', help='prefix for output file')

#----------------------CluStR----------------------------------------

def Ez(z):
    Om = 0.3
    Ov = 0.7
    return np.sqrt(Om*(1.+z)**3 + Ov)

# We'll define useful classes here
class Config:
    '''
    Used for CluStR config processing

    '''

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
    Read/Load the fits table that contains the data.

    """

    def __init__(self, cat_file_name, config):
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
    '''
    This class takes a catalog table and grabs only the relevant columns
    for the desired fit using the config dictionary.

    Config is expected to act like a dictionary
    '''

    def __init__(self, config, catalog):
        self._load_data(config, catalog)

        return

    def create_cuts(self, config, catalog):
            """
            Apply cuts to data. Will remove flags of type Boolean, Cutoff, and Range.
            """

            # Initialize an array of the same size as catalog. Elements are boolean type.
            mask = np.zeros(len(catalog), dtype=bool)

            return

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


        N = np.size(x)
        assert N == np.size(y)

        # Censored Data
        cenTF = list(config["Censored"].keys())[0]
        print(cenTF)

        if cenTF:
            cenName = config["Censored"][True]
            delta_ = catalog[cenName].astype(np.int64)

        else:
            delta_ = np.ones(N).astype(np.int64)

        if config['scale_x_by_ez'] == True:
            redshift = config['Redshift']
            x *= Ez(catalog[redshift])**(config['scaling_factor_x'])

        if config['scale_y_by_ez'] == True:
            redshift = config['Redshift']
            print(Ez(0.2)**(config['scaling_factor_y']))
            y *= (Ez(catalog[redshift]))**(config['scaling_factor_y'])


        self.x = x
        self.y = y
        self.delta_ = delta_

        return

class Fitter:
    """Runs linmix alogirthm using the regression library."""

    def __init__(self, data, config):
        """ Here we can use the super method to inherit
            the attributes from the Data class.
        """

        self.algorithm = 'linmix'
        self.data_x = data.x
        self.data_y = data.y
        self.data_xlabel = data.xlabel
        self.data_ylabel = data.ylabel
        self._constant = config['scale_line']
        self.log_data(config)
        self.fit(data)

        return

    def fit(self, data):
        '''
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma_sqr.
        '''

        self.kelly_b, self.kelly_m, self.kelly_sigsqr = reglib.run_linmix(
                                                            x=self.log_x,
                                                            y=self.log_y,
                                                            delta=data.delta_)

        self.mean_int = np.mean(self.kelly_b)
        self.mean_slope = np.mean(self.kelly_m)
        self.mean_sigsqr = np.mean(self.kelly_sigsqr)


        return

    def log_data(self, config):
        ''' Scale data to log'''

        # Log-x before pivot
        xlog = np.log(self.data_x)

        # Set pivot
        piv_type = config["piv_type"]
        if piv_type == "median":
            self.piv = np.log(np.median(self.data_x))
        else:
            self.piv = np.log(config['piv_value'])

        self.log_x = xlog - self.piv
        self.log_y = np.log(self.data_y)

        self.xmin = np.min(self.log_x)
        self.xmax = np.max(self.log_x)


        return

def main():


    #CluStR args
    args = parser.parse_args()

    config = Config(args)

    catalog = Catalog(args.cat_filename, config)

    data = Data(config, catalog)

    fitter = Fitter(data, config)

    print(f"x-pivot = {fitter.piv}")
    print(f"Mean Intercept: {np.mean(fitter.kelly_b)}")
    print(f"Mean Slope: {np.mean(fitter.kelly_m)}")
    print(f"Mean Variance: {np.mean(fitter.kelly_sigsqr)}")

    print('\n')

    print("Using Kelly Algorithm...")

    print('\nMaking Plots...')


    print('Done!')

    return

if __name__ == '__main__':
    main()
