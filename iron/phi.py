import ctypes
from itertools import chain
import sys
import os
import numpy as np
import pandas as pd
from numpy import double
from numpy.ctypeslib import ndpointer
from sklearn.model_selection import train_test_split
from iron import adjbox

from ctypes import *
phiMethods = ["extremes","range"]

'''
Generation of relevance function

@description This procedure enables the generation of a relevance function that performs a mapping between the values in a given target variable and a relevance value that is bounded by 0 (minimum relevance) and 1 (maximum relevance). This may be obtained automatically (based on the distribution of the target variable) or by the user defining the relevance values of a given set of target values - the remaining values will be interpolated.

@param y The target variable of a given data set
@param phi_parms The relevance function providing the data points where the pairs of values-relevance are known
@param method The method used to generate the relevance function (extremes or range)
@param extr_type Type of extremes to be considered: low, high or both (default)
@param control_pts Parameter required when using 'range' method, representing a 3-column matrix of y-value, corresponding relevance value (between 0 and 1), and the derivative of such relevance value
@param asym Boolean for assymetric interpolation. Default TRUE, uses adjusted boxplot. When FALSE, uses standard boxplot.


@return A dictionary with three slots with information concerning the relevance function
 \item["method"]{The method used to generate the relevance function (extremes or range)}
 \item["npts"]{?}
 \item["control_pts"]{Three sets of values identifying the target value-relevance-derivate for the first low extreme value, the median, and first high extreme value}


Example Run
df = pd.read_csv('data/accel_data.csv')
train_data, test_data = train_test_split(df, test_size=0.2, random_state=7)
ph = phi_control(train_data["acceleration"],extr_type="both")

ph = phi_control(train_data["acceleration"],extr_type="both",method = "extremes")
ph = phi_control(train_data["acceleration"],extr_type="both",method = "range")



'''

def phi_control(y, phi_parms=None, method=phiMethods,
                extr_type=None, control_pts=None, asym=True):

    if (phi_parms is not None):
        method = phi_parms["method"]
        extr_type = phi_parms["extr_type"]
        control_pts = phi_parms["control_pts"]

    if method == "range":
        control_pts = phi_range(y,control_pts=control_pts)
    else:
        # create phi_extreme and control the parameters
        control_pts = phi_extremes(y, extr_type=extr_type,asym=asym)
        method = "extremes"
    phiP = {"method": method, "npts": control_pts["npts"], "control_pts": control_pts["control_pts"]}

    return phiP

#Auxiliary function
def minmax(val_list):
    min_val = min(val_list)
    max_val = max(val_list)

    return (min_val, max_val)

'''
Relevance function for extreme target values

#@description Automatic approach to obtain a relevance function for a given target variable when the option of extremes is chosen, i.e. users are more interested in accurately predicting extreme target values

@param y The target variable of a given data set
@param extr_type Type of extremes to be considered: low, high or both (default)
@param coef Boxplot coefficient (default 1.5)
@param asym Boolean for assymetric interpolation. Default TRUE, uses adjusted boxplot. When FALSE, uses standard boxplot.



@return A dictionary with three slots with information concerning the relevance function
 \item["method"]{The method used to generate the relevance function (extremes or range)}
 \item["npts"]{?}
 \item["control_pts"]{Three sets of values identifying the target value-relevance-derivate for the first low extreme value, the median, and first high extreme value}
'''

def phi_extremes(y, extr_type="both", coef=1.5, asym=True):

    control_pts = []
    npts = None
    if asym:

        y = y.to_list()
        extr = adjbox.adjboxStats(y, coef=coef)

        r = minmax(y)
        if extr_type is None : extr_type = "both"
        if extr_type in ("both", "low"):
            ## adjL

            control_pts.append((extr["fence"][0], 1, 0))

        else:
            ## min
            control_pts.append( (r[0], 0, 0))

        ## median
        control_pts.append((extr["stats"][2], 0, 0))

        if (extr_type in ("both", "high")):

            ## adjH
            control_pts.append((extr["fence"][1], 1, 0))
        else:
            ## max
            control_pts.append( (r[1], 0, 0))

        npts = len(control_pts)

    else:

        extr = boxplot(y)

        r = minmax(y)

        if (extr_type in ("both", "low")) & (any(x<extr[0][1] for x in extr[3])):

            ## adjL
            control_pts.append( (extr[0][0], 1, 0))
        else:
            ## min
            control_pts.append( (r[0], 0, 0))

        ## median
        control_pts.append( (extr[0][2], 0, 0))

        if (extr_type in ("both", "high") )& (any(x>extr[0][4] for x in extr[3])):

            ## adjH
            control_pts.append( (extr[0][4], 1, 0))
        else:
            ## max
            control_pts.append( (r[1], 0, 0))


            npts = len(control_pts)


    latten_list = list(chain.from_iterable(control_pts))
    return {"npts": npts, "control_pts": latten_list}


'''

Custom Relevance Function

@description User-guided approach to obtain a relevance function for certain intervals of the target variable when the option of range is chosen in function phi.control, i.e. users define the relevance of values for which it is known

@param y The target variable of a given data set
@param control_pts Parameter representing a 3-column matrix of y-value, corresponding relevance value (between 0 and 1), and the derivative of such relevance value, allowing users to specify the known relevance at given target values


@return A dictionary with three slots with information concerning the relevance function
 \item["method"]{The method used to generate the relevance function (extremes or range)}
 \item["npts"]{?}
 \item["control_pts"]{Three sets of values identifying the target value-relevance-derivate for the first low extreme value, the median, and first high extreme value}



'''
def phi_range(y, control_pts) :

  if type(control_pts) is dict :
      control_pts = np.reshape(np.array(control_pts["control_pts"]),(control_pts["npts"],int(len(control_pts["control_pts"])/control_pts["npts"])))

  if (type(control_pts) is not np.ndarray) or (control_pts is None) or (np.shape(control_pts)[1] > 3) or (np.shape(control_pts)[1] < 2):
       sys.exit('The control.pts must be given as a matrix in the form: \n < x, y, m > or, alternatively, < x, y >')

  npts = len(control_pts)
  dx = control_pts[1:,0] - control_pts[0:(npts-1),0]


  if(None  in dx) or (0 in dx ) :
    sys.exit("'x' must be *strictly* increasing (non - NA)")

  if (any(x>1 for x in control_pts[:,1]))  or  any(x<0 for x in control_pts[:,1]) :
    sys.exit("phi relevance function maps values only in [0,1]")


  control_pts = control_pts[np.argsort(control_pts[:, 0])]

  if(np.shape(control_pts)[1] == 2) :
    ## based on "monoH.FC" method
    dx = control_pts[1:,0] - control_pts[0:(npts-1),0]
    dy = control_pts[1:,1] - control_pts[0:(npts-1),1]
    Sx = dy / dx
    m = (Sx[1:] + Sx[0:(npts-2),])/2
    m = np.reshape(m,(len(m),1))
    m = np.insert(m,(0,len(m)),0,axis=0)
    control_pts = np.append(control_pts,m,axis=1)


  r = minmax(y)
  npts = np.shape(control_pts)[0]
  latten_list = list(chain.from_iterable(control_pts))

  return {"npts": npts, "control_pts": latten_list}


#Auxiliary function
def phi2double(phi_parms):

    phi_parms_double = []
    if phi_parms["method"] == "extremes":
        phi_parms_double.append(0)
    elif phi_parms["method"] == "range":
        phi_parms_double.append(1)

    phi_parms_double.append(double(phi_parms["npts"]))
    phi_parms_double = np.append(phi_parms_double, phi_parms["control_pts"])

    return phi_parms_double


'''

 Obtain the relevance of data points

@description The phi function retrieves the relevance value of the values in a target variable.
It does so by resorting to the Piecewise Cubic Hermitate Interpolation Polynomial method for interpolating over 
a set of maximum and minimum relevance points. The notion of relevance is associated with rarity.
Nonetheless, this notion may depend on the domain experts knowledge

@param y The target variable of a given data set
@param phi_parms The relevance function providing the data points where the pairs of values-relevance are known
@param only_phi Boolean (default True) to return either solely the relevance values or the full data structure with the
first and second derivative the interpolated values

@return A vector or dictionary with the relevance values of a given target variable

Example Run
df = pd.read_csv('data/accel_data.csv')
train_data, test_data = train_test_split(df, test_size=0.2, random_state=7)
ph = phi_control(train_data["acceleration"],extr_type="both")
phi(test_data["acceleration"], phi_parms=ph)

Another examples for phi_control
phit = phi_control(train_data["acceleration"], extr_type="high")
phit = phi_control(train_data["acceleration"], extr_type="low")
phit = phi_control(train_data["acceleration"])


'''
def phi(y, phi_parms=None, only_phi=True):
   if phi_parms is None:
       phi_parms = phi_control(y)
   n = len(y)

   if sys.platform == "win32":
       dir = os.path.dirname(sys.modules["iron"].__file__)
       path = os.path.join(dir, "phi.dll")
       phi_c = cdll.LoadLibrary(path)
   elif  sys.platform == "darwin":
       dir = os.path.dirname(sys.modules["iron"].__file__)
       path = os.path.join(dir, "phi_mac.so")
       phi_c = cdll.LoadLibrary(path )
   elif  sys.platform == "linux":
       dir = os.path.dirname(sys.modules["iron"].__file__)
       path = os.path.join(dir, "phi_linux.so")
       phi_c = cdll.LoadLibrary(path)

   try:
       py2phi = phi_c.py2phi
       py2phi.restype = None
       py2phi.argtypes = [ctypes.c_size_t,
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

       y_phi_all = np.empty((3 * n))
       py2phi(n, y.values, phi2double(phi_parms), y_phi_all)
       phis = {"y_phi": y_phi_all[0:n], "yd_phi": y_phi_all[n:2 * n], "ydd_phi": y_phi_all[2 * n:3 * n]}

       if (only_phi):
           return phis["y_phi"]
       else:
           return phis
   except:
       print('OS %s not recognized, Only win32, macos or linux' % (sys.platform))



#Auxiliary function
def boxplot(x):
    median = np.median(x)
    upper_quartile = np.percentile(x, 75)
    lower_quartile = np.percentile(x, 25)

    iqr = upper_quartile - lower_quartile
    upper_whisker = x[x <= upper_quartile + 1.5 * iqr].max()
    lower_whisker = x[x >= lower_quartile - 1.5 * iqr].min()
    return {"lower_whisker":lower_whisker,"upper_whisker":upper_whisker,"median":median}



#Example run
def run():
    dir = os.path.dirname(sys.modules["iron"].__file__)
    df = pd.read_csv(dir + "/data/dfs_fixed.csv")
    x = df["value"]
    x = x.to_list()
    train_data, test_data = train_test_split(df, test_size=0.2, random_state=7)
    y = df["value"]
    phit = phi_control(y, extr_type="high")
    print(phit)
    phi(y, phi_parms=phit)



