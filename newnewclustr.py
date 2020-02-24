from argparse import ArgumentParser
import os
#import cPickle as pickle
from astropy.table import Table
import numpy as np
#import reglib  # Regression library
import matplotlib.pyplot as plt

class ArgumentParser:
    def __inti__(self,):

    def parse_options(self):

    parser = ArgumentParser()
    parser.add_argument('config_filename', type = str, help = 'the filename of the desired config to run')
    parser.add_argument('cat_filename', type = str, help = 'the filename of the desired catalog to open')
    #parser.add_argument('plotting_filename', type = str, help = 'the filename of the desired plotting file to run')
    valid_axes = ['1500kpc', 'lr2500', 'lr500', 'lr500cc', 't500kpc', 'tr2500', 'tr500', 'tr500cc, 'lambda']
    parser.add_argument('x', help="x-axis of plot", choices=valid_axes)
    parser.add_argument('y', help="y-axis of plot", choices=valid_axes)
    parser.add_argument('-p', help='name of output file')
    return parser.parse_arg()

class Config:

    _required_keys=[]
    _default_run_name = 'clustr'
    def __init__(self,config_file, run_options):

        with open(config_file, 'r') as stream:

            self.config_file = yaml.safe_load(stream)

        self.run_options = run_options
        if run_options.run_name is None:
            self.run_name = _default_run_name
        else:
            self.runname = run_options.run_name
        return

    def __getitem__(self,key):
        return self._config.__dict__[key]

    def __setitem(self, key, value):
        self._config._dict_[key] = value

    def __delitem__(self,key):
        del self._config.__dict__[key]

    def __contain__(slef, key):
        return key in seelf._config.__dict__

    def __len__(self):
        return len(self._config.__dict__)

    def __repr__(self):
        return repr(self._config.__dict__)

    def run_name()
pass


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

    def plot_data(self, xcol, ycol, size=8, ylog=False):
        x = self.table[xcol]
        y = self.table[ycol]

        plt.scatter(x, y)
        plt.xlabel(xcol)
        plt.ylabel(ycol)

        if ylog is True:
            plt.yscale('log')
        plt.gcf().set_size_inches(size, size) #get current figure then set size
        plt.show()

        return
    def Ez(z)
        Om = 0.3
        H_0 = 0.7
        h = H_0/100
        return np.sqrt(Om*(1.+z)**3 + h
#class Flag:
    #def __init__(self,catalog):
        #self.catalog = catalog
        #return
class Data:
    #take data frome the Catalog Class and pick rows and columns we want to fit
    #dont call flag in main call it here
    def __init__(self,config, catalog):
        self.config = config
        self.catalog = catalog
        return
    def run_config(self):
        config_results = config.rlf #run Config's function rlf and get the results. maybe?
        return config_results
    def open_catalog(self): #run Catalog's _load_catolog to have to table in this class
        table_cat = catalog._load_catolog
        return table_cat
    def get_data(self):
        #hdulist = fits.open(options.catalog)
        #data = hdulist[1].data
        #
        label_x = self.run_options[x] #how to access ?
        label_y = self.run_options[y] #
        x = self.data[label_x] #How to access?
        y = self.data[label_y] #how to access? 

        # Number of original data
        N = np.size(x)
        #where does the l come from? should it also be lambda?
        if label_x[0] == 'l' and label_y != 'lambda':
            x /= Catalog.Ez(data['redshift'])
        if label_y[0] == 'l' and label_x != 'lambda':
            y /= Catalog.Ez(data['redshift'])
ERROR STUFF

        flags = self.run_options[flags] #change. ?
        if flags is not None:
            # FIX: Should be more error handling than this!
            # FIX: Should write method to ensure all the counts are what we expect

            mask = f.create_cuts(self)
            self.x[mask] = -1
            self.y[mask] = -1

            print (
                'NOTE: `Removed` counts may be redundant, '
                'as some data fail multiple flags.'
            )

            # Take rows with good data, and all flagged data removed
            good_rows = np.all([x != -1, y != -1], axis=0)
            x = x[good_rows]
            y = y[good_rows]
            x_err = x_err[good_rows]
            y_err = y_err[good_rows]
            print 'Accepted {} data out of {}'.format(np.size(x), N)
        if np.size(x) == 0:
            print (
                '\nWARNING: No data survived flag removal. '
                'Suggest changing flag parameters in `param.config`.'
                '\n\nClosing program...\n'
            )
            raise SystemExit(2)

        print 'mean x error:', np.mean(x_err)
        print 'mean y error:', np.mean(y_err)

        hdulist.close()

        return (x, y, x_err, y_err) #put all into 1 and return 'd'

    #identify Flags
    #take flagged data out of table_cat
    #return only viable data


class Fitter:
    def __init__(self,viable_data,plotting_filename):
        self.viable_data= viable_data
        self.plotting_filename = plotting_filename
        return
    def fit(self):
        #should we use the plotting method in plotlib, write a different one in a similar file,
        #or write it directly into the code?
        pass


def main():

    args = parser.parse_args()

    config_filename = args.config_filename

    config = Config(config_filename) #(2)

    cat_file_name = args.cat_filename

    catalog = Catalog(cat_file_name, config) #(3)

    viable_data = Data(config, catalog) #(4)

    fit = Fitter.fit(viable_data) #(6)

    # Just for fun!
    catalog.plot_data('lambda', 'r500_band_lumin', ylog=True)

if __name__ == '__main__':

    main()
