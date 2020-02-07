from argparse import ArgumentParser
import os
#import cPickle as pickle
from astropy.table import Table
import numpy as np
#import reglib  # Regression library
import matplotlib.pyplot as plt
parser = ArgumentParser()
parser.add_argument(
'config_filename',
type = str,
help = 'the filename of the desired config to run'
)
parser.add_argument(
'cat_filename',
type = str,
help = 'the filename of the desired catalog to open'
)
parser.add_argument(
'plotting_filename',
type = str,
help = 'the filename of the desired plotting file to run'
)
#is it looping? what does "line" mean? what is line[0:2] pulling out?
#whats the point

class Config:
    def __init__(self,config_filename):
        self.config_filename = config_filename
        return

    def rlf(self):
        #from spencers code. Gets rid of anything commented out.
        with open(confing_filename) as config_file:
            for line in config_file:
                # Ignore empty lines and comments:
                if line[0:2] == '\n':
                    continue
                if line[0] == '#':
                    continue
                line.strip()
                line = line.split('#')[0]
                # Remove whitespace and interpret Name:Value pairs:
                line = ''.join(line.split())
                line = line.split(':')
                name, value = line[0], unicode(line[1])
        #convert strings in the param file to booleans, int, or float


class Catalog:
    #read/load the fits table that contains the data
    def __init__(self,cat_file_name,config):
        self.file_name = cat_file_name

        #self.property = config.property # for example
        #,,,

        self._load_catalog()

        return

    def _load_catalog(self):
        self.table = Table.read(self.file_name)

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
        plt.gcf().set_size_inches(size, size)
        plt.show()

        return

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
        results = config.rlf #run Config's function rlf and get the results. maybe?
        return results
    def open_catalog(self): #run Catalog's _load_catolog to have to table in this class
        table_cat = catalog._load_catolog
        return table_cat
    #identify Flags
    #take flagged data out of table_cat
    #return only viable data


class Fitter:
    def __init__(self,viable_data,plotting_filename):
        self.viable_data= viable_data
        self.plotting_filename = plotting_filename
        return
    def fit(self):

        pass


def main():

    args = parser.parse_args()
#what is "run_options" should it be an input of Config?
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
