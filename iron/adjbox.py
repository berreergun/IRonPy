import math
import os
import sys
from ctypes import cdll
import ctypes
from numpy import double
from numpy.ctypeslib import ndpointer
import numpy as np
import platform
# def fivenum(v):
#     """Returns Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for the input vector, a list or array of numbers based on 1.5 times the interquartile distance"""
#     from numpy import percentile
#     from scipy.stats import scoreatpercentile
#     try:
#         np.sum(v)
#     except TypeError:
#         print('Error: you must provide a list or array of only numbers')
#     md = np.median(v)
#     quartiles = percentile(v, [25, 50, 75])
#     return np.min(v),quartiles[0] , md,quartiles[2] , np.max(v),

def fivenum(v):
    v = np.sort(v)
    n = len(v)
    n4 = np.floor((n+3)/2) / 2
    d = [0, n4-1, (n+1)/2-1, n-n4, n-1]
    return 0.5 * (v[np.floor(d).astype(int)] + v[np.ceil(d).astype(int)])


def mcComp (x, do_reflect, do_scale, eps1, eps2, maxit = 1000, trace_lev = 1):

    if ( type(do_reflect) != bool) \
            | ( do_reflect is None) | (type(do_scale) != bool) | \
              (do_scale is None) |(not eps1 >= 0) | (not eps2 >= 0)  :
        sys.exit("Parameters are incorrect")

    ## Assumption [from caller, = mc()]: 'x' has no NAs (but can have +-Inf)
    n = int(len(x))
    eps = [eps1, eps2]
    c_iter = [maxit, trace_lev]
    if sys.platform == "win32":
        if platform.architecture()[0] == '64bit':
            dir = os.path.dirname(sys.modules["iron"].__file__)
            path = os.path.join(dir, "mc64.dll")
            mc_func = cdll.LoadLibrary(path)
        else:
            dir = os.path.dirname(sys.modules["iron"].__file__)
            path = os.path.join(dir, "mc.dll")
            mc_func = cdll.LoadLibrary(path)
    elif sys.platform == "darwin":
        dir = os.path.dirname(sys.modules["iron"].__file__)
        path = os.path.join(dir, "mc_mac.so")
        mc_func = cdll.LoadLibrary(path)
    elif sys.platform == "linux":
        dir = os.path.dirname(sys.modules["iron"].__file__)
        path = os.path.join(dir, "mc_linux.so")
        mc_func = cdll.LoadLibrary(path)


    mc_c = mc_func.mc_C
    mc_c.restype = None
    mc_c.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                     ctypes.c_size_t,
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]

    out = np.empty(1)
    iter =np.array(c_iter,dtype=np.int32)
    eps = np.array(eps)
    mc_c(np.array(x), n, eps, iter ,out, np.array(int(do_scale),dtype=np.int32))
    ans = {"medc": out[0], "eps": eps, "iter": iter, "converged": False}
    it = iter
    if it[1] == 1:
        ans["converged"] = True
    ans["iter"] = it[0]

    if do_reflect :
        out2 = np.empty(1)
        iter2 = np.array(c_iter,dtype=np.int32)
        eps2 = np.array(eps)
        mc_c(-np.array(x), n, eps2, iter2, out2, np.array(int(do_scale),dtype=np.int32))
        ans["medc2"] = out2[0]
        ans["eps2"] = eps2
        it = iter2
        if it[1] == 1:
            ans["converged2"] = True
        ans["iter2"] = it[0]

    return ans


def mc(x, na_rm = False, do_reflect = True, do_scale = True  # <- chg default to 'FALSE' ?
    , eps1 = 1e-14, eps2 = 1e-15  # << new in 0.93-2 (2018-07..)
    , maxit = 100, trace_lev = 0
    , full_result = False
    ) :

    if len(x) > 100: do_reflect = False
    ina = np.isnan(x)
    if (na_rm):
        x = [y for y in x if np.isnan(y) == False]
    elif (False in ina) :
        return ("NA_real_")
    ## ==> x is NA-free from here on

    ## if(length(l.. <- list(...)))
    ##     stop("In mc(): invalid argument(s) : ",
    ##          paste(sQuote(names(l..)), collapse=","), call. = FALSE)
    rr = mcComp(x, do_reflect, do_scale=do_scale, eps1=eps1, eps2=eps2,
                  maxit=maxit, trace_lev = trace_lev)

    if (not (rr["converged"])) | (do_reflect and not (rr["converged2"])):
        sys.exit("mc(): not 'converged'")

    if (do_reflect) :
        m = (rr["medc"] - rr["medc2"]) / 2
    else :
        m = rr["medc"]

    if (full_result):
       return (m,rr)
    else :
       return m


def minmax(val_list):
    min_val = min(val_list)
    max_val = max(val_list)
    return (min_val, max_val)

def adjboxStats (x, coef = 1.5, a = -4, b = 3,do_conf = True, do_out = True):

    if coef < 0 :
       sys.exit("'coef' must not be negative")
    nna = [True for y in x if y is not None]
    n = sum(nna)# including +/- Inf
    stats = np.array(fivenum(x))
    iqr = stats[3] - stats[1]
    fence = [None,None]
    if coef == 0:
          do_out = False # no whiskers to be drawn
    else : ## coef > 0
        if iqr is not None:
            medc = mc(x, na_rm=True)
            if medc >= 0:
                fence = [stats[1] - (coef * math.exp(a * medc) * iqr), stats[3] + coef * math.exp(b * medc) * iqr]
            else:
                fence = [stats[1] - coef * math.exp(-b * medc) * iqr,stats[3] + coef * math.exp(-a * medc) * iqr]

            out =  (x < fence[0]) | (fence[1] < x )

        else : out =  math.isfinite(x)
        if any(out[nna]):
            stats[0] = minmax(np.array(x)[np.invert(out)])[0]
            stats[4] = minmax(np.array(x)[np.invert(out)])[1]

    if (do_conf):
            conf = stats[2] + np.array((-1.58,1.58)) * (iqr/math.sqrt(n))
    if (do_out):
        res = {"stats" : stats, "n" : n, "conf" : conf, "fence" : fence,"out" : np.array(x)[out & nna]}
    else :
        res = {"stats" : stats, "n" : n, "conf" : conf, "fence" : fence,"out" : 0}


    return res