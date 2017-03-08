import linmix # Kelly algorithm package ported to Python
import numpy as np
import numpy.random as npr
from scipy import stats
import scipy.optimize as sop
from inputParameters import beta1, beta2

npr.seed(800)
def scatter_cal(x,y,slope,intercept,dof):
   sig2 = sum((np.array(y) - (slope*np.array(x)+intercept))**2) / dof
   return np.sqrt(sig2)


def invScalingRelation(tInt,tSlope,tSig):
   xs = 1.0 / (1.0 + beta2*(tSig**2)/(tSlope**2))
   invInt = xs * ( - tInt / tSlope + beta1*(tSig**2)/(tSlope**2) )
   invSlope = xs / tSlope
   invSig = np.sqrt(xs * (tSig**2) / (tSlope**2) )
   return invInt, invSlope, invSig


def ninvScalingRelation(tInt,tSlope,tSig):
   invInt = ( - tInt / tSlope )
   invSlope = 1.0 / tSlope
   invSig = np.sqrt( (tSig**2) / (tSlope**2) )
   return invInt, invSlope, invSig


def obsScalingRelation(tInt1,tSlope1,tSig1,tInt2,tSlope2,tSig2,r):
   # First order approximation
   invInt1 = ( - tInt1 / tSlope1 + beta1*(tSig1**2)/(tSlope1**2) )
   invSlope1 = 1.0 / tSlope1
   invSig1 = np.sqrt( (tSig1**2) / (tSlope1**2) )
   invSig2 = np.sqrt( (tSig2**2) / (tSlope2**2) )

   x1 = 1.0 / (1.0 + beta2*invSig1**2)
   inter = tInt2 + x1*tSlope2*( invInt1 \
                                - (r * invSig1 * invSig2) \
                                   * ( beta1 + beta2 * tInt1 / tSlope1) )
   slope = x1 * tSlope2 * ( invSlope1 \
                            + beta2 * r * invSig1 * invSig2 / tSlope1 )
   sig = tSlope2 * np.sqrt(x1) *\
           np.sqrt( invSig2**2 + invSig1**2 - 2*r*invSig1*invSig2\
                    + beta2*invSig1**2*invSig2**2*(1.-r**2) )
   return inter, slope, sig


def nobsScalingRelation(tInt1,tSlope1,tSig1,tInt2,tSlope2,tSig2,r):
   # First order approximation
   invInt1 = ( - tInt1 / tSlope1 )
   invSlope1 = 1.0 / tSlope1
   invSig1 = np.sqrt( (tSig1**2) / (tSlope1**2) )
   invSig2 = np.sqrt( (tSig2**2) / (tSlope2**2) )

   inter = tInt2 + tSlope2*( invInt1 )
   slope = tSlope2 * ( invSlope1 )
   sig = tSlope2 * np.sqrt( invSig2**2 + invSig1**2 - 2*r*invSig1*invSig2 )
   return inter, slope, sig

def findY(Y,invSig):
   xs = 1.0 / (1.0 + beta2*Y**2)
   f = invSig - np.sqrt(xs * Y**2 )
   return f

def solveForZ_old(Z,Y,sigZY,slopeZY,ySlope,r):
   xsy = 1.0 / (1.0 + beta2*Y**2)
   slopeZ = slopeZY * ySlope / xsy / (1.0 + r*beta2*Y*Z)
   f = sigZY**2 - slopeZ**2 * xsy * \
        ( Y**2 + Z**2 - 2.*r*Y*Z + beta2*(Y**2)*(Z**2)*(1.-r**2) )
   return f

def solveForZ(Y,sigZY,slopeZY,ySlope,r):
   p0 = slopeZY**2*ySlope**2*(1.0 + beta2*Y**2*(1.-r**2))
   p1 = -slopeZY**2*ySlope**2*2.*r*Y - sigZY**2*beta2*r*Y
   p2 = slopeZY**2*ySlope**2*Y**2 - sigZY**2
   Z1,Z2 = np.roots([p0,p1,p2])
   if np.iscomplex(Z1): return 0.,0.
   return Z1,Z2

# calculate the true intercept, slope, and scatter of inverse of scaling
# relation assuming beta1 and beta2 is known (E14 notation)
def inferScalingRelationThroughInverse(infInt,infSlope,infSig):
   Y = sop.fsolve(findY,infInt/infSlope,args=infSig)[0] #sig / slope
   xs = 1.0 / (1.0 + beta2*Y**2)
   Slope = xs / infSlope
   Scatter = Y * Slope
   Intercept = - Slope * (infInt / xs - beta1 * Y**2)
   return Intercept, Slope, Scatter #OK


def inferScalingRelationThroughHidenVaribale(\
         infInt, infSlope, infSig, yInt, ySlope, ySig, r, gInt, gSlope, gSig,\
         Zg=0.0):
   Y = ySig / ySlope #sig / slope
   xsy = 1.0 / (1.0 + beta2*Y**2)
   #Z = gSig / gSlope #initial guess
   #Z = sop.fsolve(solveForZ,Z,args=(Y,infSig,infSlope,ySlope,r))[0]
   Z1,Z2 = solveForZ(Y,infSig,infSlope,ySlope,r)
   if (Z1 > Z2 ): Z = Z1
   else: Z = Z2
   #if (Zg <= 0.0): Z = Z1
   #else: Z = Z2
   #if (Z1 <= 0.0):
   #   if (Z2 <= 0.0): Z = 0.
   #   else: Z = Z2
   #else:
   #   if (Z2 <= 0.0): Z = Z1
   #   else:
   #      if (Z1 > Z2): Z = Z1
   #      else: Z = Z2
   Slope = infSlope * ySlope / xsy / (1.0 + r*beta2*Y*Z)
   Scatter = Z * Slope
   invyInt = ( - yInt/ySlope + beta1*(ySig**2)/(ySlope**2) )
   Intercept = infInt - Slope*xsy*(invyInt - r*Y*Z*(beta1 + beta2*yInt/ySlope))
   return Intercept, Slope, Scatter, Z #OK

   #Y = ySig / ySlope #sig / slope
   #xsy = 1.0 / (1.0 + beta2*Y**2)
   #Z = sop.fsolve(solveForZ,-10.0,args=(Y,infSig,infSlope,ySlope,r))[0]
   #Slope1 = infSlope * ySlope / xsy / (1.0 + r*beta2*Y*Z)
   #Scatter1 = Z * Slope
   #Z = sop.fsolve(solveForZ,5.,args=(Y,infSig,infSlope,ySlope,r))[0]
   #Slope = infSlope * ySlope / xsy / (1.0 + r*beta2*Y*Z)
   #Scatter = Z * Slope
   #invyInt = ( - yInt/ySlope + beta1*(ySig**2)/(ySlope**2) )
   #Intercept = infInt - Slope*xsy*(invyInt - r*Y*Z*(beta1 + beta2*yInt/ySlope))
   #return Intercept, Slope, Scatter #OK



#def makeLinearRegression(xObs,yObs,xerr,yerr):
#   print len(xObs), len(yObs), len(xerr), len(yerr)
#   delta = np.ones(len(xerr)); xycov = np.zeros(len(xerr))
#   model = linmix.LinMix(xObs,yObs,xerr,yerr,xycov,delta,2,2)


"""
    Args:
        x(array_like): The observed independent variable.
        y(array_like): The observed dependent variable.
        xsig(array_like): 1-sigma measurement errors in x.
        ysig(array_like): 1-sigma measurement errors in y.
        xycov(array_like): Covariance between the measurement errors in x and y.
        delta(array_like): Array indicating whether a data point is
                           censored (i.e., not detected), or not.
                           If delta[i] == 1, then the ith source is detected.
                           If delta[i] == 0, then the ith source is not
                           detected and y[i] will be interpreted as an upper
                           limit.  Note that if there are censored data points,
                           then the maximum-likelihood estimate (alpha, beta,
                           sigsqr) is not valid. By default,
                           all data points are assumed to be detected.
        K(int): The number of Gaussians to use in the mixture model
                for the distribution of xi.
        nchains(int): The number of Monte Carlo Markov Chains to instantiate.
"""

def makeLinearRegression(xObs,yObs,xerr,yerr):
   print len(xObs), len(yObs), len(xerr), len(yerr)
   delta = np.ones(len(xerr)); xycov = np.zeros(len(xerr))
   model = linmix.LinMix(xObs,yObs,xerr,yerr,xycov,delta,2,2)
   model.run_mcmc(5000, 10000, silent=False)
   # return intercept, slope, scatter
   return model.chain['alpha'], model.chain['beta'],\
          np.sqrt(model.chain['sigsqr'])

def makeOLR(x,y):
   slope, intercept, r_value, p_value, _ = stats.linregress(x,y)
   sig = scatter_cal(x,y,slope,intercept,len(x)-2)
   return intercept, slope, sig
