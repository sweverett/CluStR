import argparse as args
import os
#import cPickle as pickle
import astropy.io.fits as fits
import numpy as np
#import reglib  # Regression library
#import plotlib

class Config:
    def __init__(self,config_filename):
        pass
class Catalog:
    #read/load the fits table that contains the data
    def __init__(self,fits_file):
        self.fits = fits_file
        pass
class Flag:
    def __init__(self,data):
        self.data = data
        pass
class Data:
    #take data frome the Catalog Class and pick rows and columns we want to fit
    def __init__(self,data,flags):
        self.data = data
        self.flags = flags
        pass
class Fitter:
    def __init__(self,data_set):
        self.data_set = data_set
    def fit(self):
        pass


def main():

#what is "run_options" should it be an input of Config?
    config_filename = 'param.config'

    config_file = Config(config_filename) #(2)

    file_name = fits.open("y3a2-6.4.22+2_peak_fixed_maybe_final_merged(2).fits")

    full_data_set = Catalog(file_name) #(3)

    flagged_data = Flag(full_data_set) #(4)

    viable_data = Data(full_data_set,flagged_data) #(4)

    Fitter.fit(viable_data) #(6)




if __name__ == '__main__':

    main()
