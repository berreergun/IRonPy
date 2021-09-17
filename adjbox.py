import math
import sys
from ctypes import cdll
import ctypes
from numpy import double
from numpy.ctypeslib import ndpointer
import numpy as np
def fivenum(v):
    """Returns Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for the input vector, a list or array of numbers based on 1.5 times the interquartile distance"""

    from scipy.stats import scoreatpercentile
    try:
        np.sum(v)
    except TypeError:
        print('Error: you must provide a list or array of only numbers')
    q1 = scoreatpercentile(v,25)
    q3 = scoreatpercentile(v,75)
    iqd = q3-q1
    md = np.median(v)
    whisker = 1.5*iqd
    return np.min(v), md-whisker, md, md+whisker, np.max(v),



def mcComp (x, do_reflect, do_scale, eps1, eps2, maxit = 1000, trace_lev = 1):

    '''stopifnot(is.logical(do_reflect), length(do_reflect) == 1L, !is.na(do_reflect),
              is.logical(do_scale),   length(do_scale)   == 1L, !is.na(do_scale),
              is.1num(eps1), eps1 >= 0,
              is.1num(eps2), eps2 >= 0,
              length(maxit     <- as.integer(maxit)) == 1,
              length(trace.lev <- as.integer(trace.lev)) == 1
              )'''

    ## Assumption [from caller, = mc()]: 'x' has no NAs (but can have +-Inf)
   # x <- as.numeric(x)
    n = int(len(x))
    # todo double olmalı
    eps = ([eps1, eps2])
    c_iter = (maxit, trace_lev)
    ## NAOK=TRUE: to allow  +/- Inf to be passed

    phi_c = cdll.LoadLibrary(r"C:\Users\MONSTER\PycharmProjects\IRonProject\src\mc.dll")
    py2phi = phi_c.py2phi
    py2phi.restype = None
    py2phi.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                       ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]

    out = double(1)
    py2phi(x, n, eps,c_iter,out,do_scale)


def mc(x, na_rm = False, do_reflect = True, do_scale = True  # <- chg default to 'FALSE' ?
    , eps1 = 1e-14, eps2 = 1e-15  # << new in 0.93-2 (2018-07..)
    , maxit = 100, trace_lev = 0
    , full_result = False
    ) :

    if len(x) > 100: do_reflect = False

    #x = as.numeric(x)
    ina = x.isna()
    if (na_rm):
        x = x[not ina]
    elif (False in ina) :
        return ("NA_real_")
    ## ==> x is NA-free from here on

    ## if(length(l.. <- list(...)))
    ##     stop("In mc(): invalid argument(s) : ",
    ##          paste(sQuote(names(l..)), collapse=","), call. = FALSE)
    rr = mcComp(x, do_reflect, do_scale=do_scale, eps1=eps1, eps2=eps2,
                  maxit=maxit, trace_lev = trace_lev)

    if (not (rr[["converged"]])) | (do_reflect and not(rr[["converged2"]])) :
        sys.exit("mc(): not 'converged'")

    '''
          stop("mc(): not 'converged' ",
          if (!conv1) paste("in", rr[["iter"]], "iterations"),
          if (do_reflect & & !conv2)
          paste( if (!conv1)" *and*",
          "'reflect part' in", rr[["iter2"]], "iterations"),
          "; try enlarging eps1, eps2 !?\n")
          '''


    if (do_reflect) :
        m = (rr[["medc"]] - rr[["medc2"]]) / 2
    else :
        m = rr[["medc"]]

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
    nna = not (None in x)
    n = sum(nna)# including +/- Inf
    stats = fivenum(x)
    iqr = stats[4] - stats[2]
    fence = [None,None]
    if coef == 0:
          do_out = False # no whiskers to be drawn
    else : ## coef > 0
        if (None in iqr)  :
            medc = mc(x,na_rm = True)
            if (medc >= 0) :
                  fence = [stats[2] - coef * a ** medc * iqr,stats[4] + coef * b ** medc * iqr]
            else :
                  fence = [stats[2] - coef * -b ** medc * iqr,stats[4] + coef * -a ** medc * iqr]
            out = x < fence[1] | fence[2] < x

        else : out =  math.isfinite(x)
        #todo any(out[nna], na_rm = True) ne tam?
        if any(out[nna], na_rm = True):
            #todo minmax na.rm parametresi mne için
            stats[1] = minmax(x[not out])[0]
            stats[5] = minmax(x[not out])[5]
              #stats[c(1, 5)] <- minmax(x[not out])

    if (do_conf):
        #todo nedir?
            conf = stats[3] + [-1.58, 1.58] * iqr/math.sqrt(n)
    #list(stats = stats, n = n, conf = conf, fence = fence,
    if (do_out):
        res = {"stats" : stats, "n" : n, "conf" : conf, "fence" : fence,"out" : x[out and nna]}
    else :
        res = {"stats" : stats, "n" : n, "conf" : conf, "fence" : fence,"out" : 0}


    return res