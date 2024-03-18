import numpy as np
import pandas as pd
from iron import phi

def sera_loss(ph, step=.001):
    def sera_loss_inner(trues: np.array, preds: np.array):
        trues = pd.Series(trues, dtype='float64')
        preds = pd.Series(preds, dtype='float64')
        phis = np.array(phi.phi(trues, ph), dtype=np.float64)
        grad = sera_first_deriv(trues, preds, phis, step)
        hess = sera_second_deriv(trues, preds, phis, step)
        return grad, hess
    return sera_loss_inner

def sera_first_deriv(trues: pd.Series, preds: pd.Series, phis, step):
    if trues.dtype != 'float64':
        trues = trues.astype('float64')
    if preds.dtype != 'float64':
        preds = preds.astype('float64')
    th = np.arange(0,1+step,step)
    errors = np.zeros((len(trues), len(th)))
    for i in range(len(th)):
        errors[:, i] = first_deriv_helper(trues, preds, phis, th[i], True)
    areas = np.zeros((len(th)-1, len(trues)))
    for j in range(len(th) - 1):
        areas[j, :] = step * (errors[:, j] + errors[:, j + 1]) / 2
    return 2 * areas.sum(axis=0)

def first_deriv_helper(trues, preds, phi_trues, th, include_equals=True):
    if include_equals:
        zr = np.ones(trues.shape[0])
    else:
        zr = np.zeros(trues.shape[0])
    mask = phi_trues - th
    mask = np.heaviside(mask, zr)    
    diff = (preds - trues) * mask
    return diff

def sera_second_deriv(trues: pd.Series, preds: pd.Series, phis, step):
    if trues.dtype != 'float64':
        trues = trues.astype('float64')
    if preds.dtype != 'float64':
        preds = preds.astype('float64')
    th = np.arange(0,1+step,step)
    errors = np.zeros((len(trues), len(th)))
    for i in range(len(th)):
        errors[:, i] = second_deriv_helper(trues, phis, th[i], True)
    areas = np.zeros((len(th)-1, len(trues)))
    for j in range(len(th) - 1):
        areas[j, :] = step * (errors[:, j] + errors[:, j + 1]) / 2
    return 2 * areas.sum(axis=0)

def second_deriv_helper(trues, phi_trues, th, include_equals=True):
    if include_equals:
        zr = np.ones(trues.shape[0])
    else:
        zr = np.zeros(trues.shape[0])
    mask = phi_trues - th
    mask = np.heaviside(mask, zr)
    return mask