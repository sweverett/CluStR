from argparse import ArgumentParser
import os
import pickle as cPickle
from astropy.table import Table
import numpy as np
import reglib  # Regression library
import matplotlib.pyplot as plt
import linmix
import yaml

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
        self.config_filename = args.config_file
        self.x = args.x
        self.y = args.y

        with open(self.config_filename, 'r') as stream:
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

    #def run_name(self):
    #    pass

class Catalog:
    #read/load the fits table that contains the data
    def __init__(self,cat_file_name,config):
        self.file_name = cat_file_name

        #self.property = config.property # for example

        self._load_catalog()

        return

    def _load_catalog(self):
        self.cat_table = Table.read(self.file_name)

        # could do other things...

        return


def Ez(z):
    Om = 0.3
    H_0 = 0.7
    h = H_0/100
    return np.sqrt(Om*(1.+z)**3 + h)

class Data:
    '''
    This class takes a catalog table, and grabs only the relevant columns
    for the desired fit using the config dictionary.

    config is expected to act like a dictionary
    '''
    #take data frome the Catalog Class and pick rows and columns we want to fit
    #dont call flag in main call it here
    def __init__(self, config, catalog):
        self.get_data(config, catalog)

        return

    def get_data(self, config, catalog):
        '''
        Obtains x, y, x errors, and y errors from config & catalog files. 
        '''
        self.xlabel = fits_label(config.x)
        self.ylabel = fits_label(config.y)
        x = catalog.cat_table[self.xlabel]
        y = catalog.cat_table[self.ylabel]

        # Number of original data
        #N = np.size(x)

        # Scale data if a luminosity
        if config['scale_x_by_ez']:
            x /= Ez(catalog.cat_table['redshift'])
        if config['scale_y_by_ez']:
            y /= Ez(catalog.cat_table['redshift'])

        self.x_err = (catalog.cat_table[self.xlabel+'_err_low'] + catalog.cat_table[self.xlabel+'_err_high']) / 2.
        self.y_err = (catalog.cat_table[self.ylabel+'_err_low'] + catalog.cat_table[self.ylabel+'_err_high']) / 2.

        # Take rows with good data, and all flagged data removed
        good_rows = np.all([x != -1, y != -1], axis=0)

        self.x = x[good_rows]
        self.y = y[good_rows]
        self.x_err = self.x_err[good_rows]
        self.y_err = self.y_err[good_rows]

        #print(self.xlabel)
        #print(self.ylabel)
        #print('mean x error:', np.mean(self.x_err))
        #print('mean y error:', np.mean(self.y_err))

        return (self.xlabel, self.ylabel, self.x, self.y, self.x_err, self.y_err)

    # Plotting the x and y data.
    def plot_data(self, x_axis, y_axis):
        '''
        Plots x-axis and y-axis data.
        '''
        plt.scatter(x_axis, y_axis)
        plt.title('{} vs. {}'.format(self.xlabel, self.ylabel))
        plt.xlabel('{}'.format(self.xlabel))
        plt.ylabel('{}'.format(self.ylabel))
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.grid() #add (True, which='both') to display major and minor grids
        plt.show()

        return



class Fitter:
    def __init__(self, viable_data):
        self.fit(viable_data)
        #self.viable_data = data
        #self.plotting_filename = plotting_filename

        return

    def fit(self, viable_data):
        '''
        Calculates fit parameters using the Kelly method (linmix) and returns
        intercept, slope, and sigma.
        '''
        x_obs = viable_data[2]
        y_obs = viable_data[3]
        x_err = viable_data[4]
        y_err = viable_data[5]

        #run linmix
        print("Using Kelly Algorithm...")
        kelly_b, kelly_m, kelly_sig = reglib.run_linmix(x_obs, y_obs, x_err, y_err)

        return(x_obs, y_obs, x_err, y_err), (kelly_b, kelly_m, kelly_sig)

    def scale(self, x_obs, y_obs, x_err, y_err):
        ''' Scale data for plotting'''
        log_x = np.log(x_obs)
        x_piv = np.median(log_x)
        log_y = np.log(y_obs)

        return [log_x-x_piv, log_y, x_err/x_obs, y_err/y_obs, x_piv]

"""class SaveData(Fitter):
    def __init__(self, run_options, parameters)
"""

def main():

    args = parser.parse_args()

    config = Config(args) #(2)

    cat_file_name = args.cat_filename

    catalog = Catalog(cat_file_name, config) #(3)

    data = Data(config, catalog) #(4)

    viable_data = data.get_data(config, catalog) #check that this is the correct way to access

    #what would to plotting filename be i put a placepholder
    #plot_filename = args.plotting_filename

    #fitter = Fitter(viable_data) #(6)

    #linmixfit = fitter.fit(viable_data)
    
    # Scatter plot
    x_axis = data.x
    y_axis = data.y
    data.plot_data(x_axis,y_axis)

if __name__ == '__main__':
    main()
