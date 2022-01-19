"""Luminosity Library for CluStR"""

from astropy.table import Table
import numpy as np

from clustr import Data

data1 = Table.read('Catalog/y3a2-6.4.22+2_v3.2.fits', format='fits')

def detectedWithTemp(data, x, y, detected='Detected', temp='r2500_temperature'):
  """ 
  This function will return observations 
  that are Detected and have temperature values.
  
  Will return x and y values. 
  """
  data[temp] = np.ma.filled(data[temp], fill_value=0)

  firstCond = data[temp] != 0
  secondCond = data[detected] == True
  yWithValues = y[firstCond & secondCond]
  xWithValues = x[firstCond & secondCond]

  return xWithValues, yWithValues

def detectedWithNoTemp(data, x, y, xerr, yerr, detected='Detected', temp='r2500_temperature'):
  """ 
  This function will return observations 
  that are Detected and have NO Temperature values.
  
  Will return x and y values. 
  """

  # Some values are empty. We'll use the filled() method to set them to zero.
  data1[temp] = data1[temp].filled(0)

  firstCond = data1[temp] == 0
  secondCond = data1[detected] == True

  yWithValues = y[firstCond & secondCond]
  yErrWithValues = yerr[firstCond & secondCond]
  xWithValues = x[firstCond & secondCond]
  xErrWithValues = xerr[firstCond & secondCond]

  

  return xWithValues, yWithValues, xErrWithValues, yErrWithValues

def uppLimUndetected(data, x, y, detected='Detected'):
  """ 
  This function will return observations 
  that are Not Detected and have upper-limits on luminosity values.
  
  Will return x and y values. 
  """

  secondCond = data[detected] == False
  yWithValues = y[secondCond]
  xWithValues = x[secondCond]

  return xWithValues, yWithValues

