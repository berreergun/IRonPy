import numpy as np
from iron import phi
import pandas as pd
from matplotlib import pyplot as plt

def sera(trues, preds, phi_trues=None, ph=None, pl=False, step = 0.001, return_err = False, include_equals=True):
    if ((phi_trues is None) & (ph is None)):
        raise ValueError("You need to input either the parameter phi_trues or ph.")
    if phi_trues is None :
       phi_trues = autophi(trues, ph)

    if preds.ndim > 1:
        if preds.shape[1] > 1:
            return __sera_multiple(trues, preds, phi_trues=phi_trues, ph=ph, pl=pl, step=step, return_err=return_err, include_equals=include_equals)

    trues = trues.values

    sera_v = sera_t_vectorizer(trues, preds, phi_trues, include_equals)
    #generate a list of points for sert eval
    th = np.arange(0, 1+step, step)
    #eval sert on candidate t's
    errors = sera_v(th)
    #use trapezoid to numerically compute integral
    res = np.trapz(errors, th)
    if pl:
        plt.plot(th, errors)
        plt.ylabel("Error")
        plt.xlabel("Threshold")
        plt.show()

    if return_err :
       return {"sera": res, "errors": errors, "thrs": th}
    else:
       return res

def sera_t_vectorizer(trues, preds, phi_trues, include_equals=True):   
    def local_sera(th):
        return sera_helper(trues, preds, phi_trues, th, include_equals)
    return np.vectorize(local_sera)

def sera_helper(trues, preds, phi_trues, th, include_equals=True):
    if include_equals:
        zr = np.ones(trues.shape[0])
    else:
        zr = np.zeros(trues.shape[0])
    mask = phi_trues - th
    mask = np.heaviside(mask, zr)
    sq_err = (preds - trues) ** 2
    sq_err = sq_err * mask
    return np.sum(sq_err)

def autophi(y, ph, out_dtype=np.float64):
    y = np.float64(np.array(y))
    phi_values = phi.phi(pd.Series(y), ph)
    return np.array(phi_values, dtype=out_dtype)

def __sera_multiple(trues, preds, phi_trues=None, ph=None, pl=False, step = 0.001, return_err = False, include_equals=True):
    res = {}
    thrs = np.arange(0, 1+step, step)
    errors = np.zeros((len(thrs), preds.shape[1]), dtype=float)
    i = 0
    for col in preds:
        temp_sera = sera(trues, preds.loc[: , col], phi_trues=phi_trues, ph=ph, step=step, return_err=True, include_equals=include_equals)
        res[col] = temp_sera['sera']
        errors[:, i] = temp_sera['errors']
        i += 1
        if pl:
            plt.plot(thrs, temp_sera['errors'], label=col)
    res = pd.DataFrame.from_dict(res, orient='index')[0]
    if pl:
        plt.ylabel("Error")
        plt.xlabel("Threshold")
        plt.legend()
        plt.show()        
    if return_err:
        return {"sera": res, "errors": errors, "thrs": thrs}
    return res