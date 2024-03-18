import numpy as np

from iron import phi
import pandas as pd
from matplotlib import pyplot as plt

def sera (trues, preds, phi_trues = None, ph = None, pl = False, m_name = "Model", step = 0.001,return_err = False) :
    if not isinstance(preds, pd.DataFrame):
        preds = pd.DataFrame(preds)
    if ((phi_trues is None) & (ph is None)):
        raise ValueError("You need to input either the parameter phi_trues or ph.")
    if phi_trues is None :
       phi_trues = phi.phi(trues, ph)
    trues = trues.values
    tbl = pd.DataFrame(
        {'trues': trues,
         'phi_trues': phi_trues,
         })
    tbl = pd.concat([tbl, preds], axis=1)
    ms = list(tbl.columns[2:])
    th = np.arange(0, 1 + step, step)
    errors = []
    for ind in th:
        errors.append(
            [sum(tbl.apply(lambda x: ((x['trues'] - x[y]) ** 2) if x['phi_trues'] >= ind else 0, axis=1)) for y in ms])

    areas = []
    for x in range(1, len(th)):
        areas.append([step *(errors[x - 1][y] + errors[x][y]) / 2 for y in range(len(ms))])
    areas = pd.DataFrame(data=areas, columns=ms)
    res = areas.apply(lambda x: sum(x))
    names = res.index
    if pl:
        plt.plot(th, errors, label=names)
        plt.ylabel("Error")
        plt.xlabel("Threshold")
        if len(names) > 1:
            plt.legend()
        plt.show()

    if return_err :
       return {"sera":res, "errors":errors, "thrs" :th}
    else:
       return res