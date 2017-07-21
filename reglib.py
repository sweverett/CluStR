''' Regression library for CluStR '''

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy as np
import linmix

# Imports the necessary R packages needed to run lrgs in python
RLRGS = importr('lrgs')  # Multivariate regression package by Adam Mantz

# Set some aliases for useful R functions
RARRAY = robjects.r('array')
RMATRIX = robjects.r('matrix')
RNORM = robjects.r('rnorm')
RC = robjects.r('c')
RLM = robjects.r('lm')

# pylint: disable=invalid-name

def run_lrgs(x, y, err_x, err_y, _xycov=None, nmc=500, dirichlet=True):
    '''
    Runs the lrgs regression algorithm written in R by interfacing through
    rpy2. For our purposes, inputs should be in scaled (log) form. (For the
    moment, only works for on-diagonal elements of the covariance matrix.) nmc
    is the length of the markov chain.
    '''

    # pylint: disable = too-many-arguments
    # pylint: disable = too-many-locals

    # Make sure dimensions are correct
    assert np.size(x) == np.size(y)
    assert np.size(err_x) == np.size(err_y)
    assert np.size(x) == np.size(err_x)

    # Convert x and y to r vectors
    rx = robjects.FloatVector(x)
    ry = robjects.FloatVector(y)
    rx_err = robjects.FloatVector(err_x)
    ry_err = robjects.FloatVector(err_y)

    # Set up covariance matrix
    M = RARRAY(0.0, dim=RC(2, 2, np.size(rx)))

    for i in range(np.size(rx)):
        M.rx[1, 1, i+1] = rx_err[i]
        M.rx[2, 2, i+1] = ry_err[i]

    # Set some R equivalents
    TRUE = robjects.BoolVector([True])
    FALSE = robjects.BoolVector([False])

    if dirichlet:
        d = TRUE
    else:
        d = FALSE

    # Run MCMC
    posterior = RLRGS.Gibbs_regression(rx, ry, M, nmc, dirichlet=d,
                                       trace='bsg', mention_every=50)

    # Extract relevant data from posterior
    B = np.array(posterior[0])  # Parameter chain
    S = np.array(posterior[1])[0][0]
    # ^ Scatter chain (only intrinsic scatter for the moment!)

    # Prepare lrgs fit chains
    intercept = B[0][0]
    slope = B[1][0]
    sigma = np.sqrt(S)

    # Return fit parameters consistently with run_linmix
    return (intercept, slope, sigma)


def run_linmix(x, y, err_x, err_y, Nmin=5000, Nmax=10000, vb=True):
    # pylint: disable = too-many-arguments
    ''' Runs the Kelly regression algorithm through the package linmix.'''

    ''' For convenience, here are the linmix arguments:

        Linmix Args:
            x(array_like): The observed independent variable.
            y(array_like): The observed dependent variable.
            xsig(array_like): 1-sigma measurement errors in x.
            ysig(array_like): 1-sigma measurement errors in y.
            xycov(array_like): Covariance between the measurement errors in x
                               and y.
            delta(array_like): Array indicating whether a data point is
                               censored (i.e., not detected), or not.
                               If delta[i] == 1, then the ith source is
                               detected. If delta[i] == 0, then the ith source
                               is not detected and y[i] will be interpreted as
                               an upper limit. Note that if there are censored
                               data points, then the maximum-likelihood
                               estimate (alpha, beta, sigsqr) is not valid. By
                               default, all data points are assumed to be
                               detected.
            K(int): The number of Gaussians to use in the mixture model
                    for the distribution of xi.
            nchains(int): The number of Monte Carlo Markov Chains to
                          instantiate.
    '''


    # Make sure dimensions are correct
    assert np.size(x) == np.size(y)
    assert np.size(err_x) == np.size(err_y)
    assert np.size(x) == np.size(err_x)

    L = np.size(x)

    # FIX: Implement censored data!
    # Run linmix MCMC
    delta = np.ones(L)
    xycov = np.zeros(L)
    model = linmix.LinMix(x, y, err_x, err_y, xycov, delta, 2, 2)
    model.run_mcmc(Nmin, Nmax, silent=vb)

    # return intercept, slope, intrinsic scatter
    intercept = model.chain['alpha']
    slope = model.chain['beta']
    sigma = np.sqrt(model.chain['sigsqr'])

    # Return fit parameters consistently with run_lrgs
    return (intercept, slope, sigma)


def check_convergence(_intercept, _slope, _sigma):
    '''
    FIX: In future, should implement a function that checks the convergence of
    the MCMC using autocorrelation, etc. and display plots/statistics.
    '''
    pass
