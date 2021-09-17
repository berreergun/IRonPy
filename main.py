# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import ctypes
from ctypes import *
import pandas as pd
import numpy as np
from numpy.ctypeslib import ndpointer

import sera
from sklearn.cross_validation import train_test_split
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.25, random_state = 0)



trues = (2,3,36,63636,669,6)
phi_trues= (12,13,136,163636,1669,16)
preds1 = (12,13,136,163636,1669,16)
preds2 = (23,33,363,636363,6639,63)

preds = pd.DataFrame({"preds1":preds1,"preds2":preds2})
preds = (23,33,363,636363,6639,63)
print(preds)
step = 0.001

res =sera.sera(preds=preds,phi_trues=phi_trues,trues=trues)

print(res)

tbl = pd.DataFrame(
    {'trues': trues,
     'phi_trues': phi_trues,
    })
tbl = pd.concat([tbl, preds], axis=1)

# Creates DataFrame.

ms = list(tbl.columns[2:])

print(ms)
#preds_size = len(tbl['preds'][0])
th = np.arange(0, 1 + step, step)
errors = []

# todo ne biçim saturdur
for ind in th:
       errors.append([sum(tbl.apply(lambda x: ((x['trues'] - x[y]) ** 2) if x['phi_trues'] >= ind else 0, axis=1)) for y in ms])

areas = []
for x in range(1, len(th)):
    areas.append([(errors[x - 1][y] + errors[x][y]) / 2 for y in range(len(ms))])

areas = pd.DataFrame(data = areas, columns=ms)
print(areas)
res = areas.apply(lambda x: sum(x))
print(res)

step = 0.001
th = np.arange(0, 1+ step, step)
print(th)

phi_c = cdll.LoadLibrary(r"C:\Users\MONSTER\PycharmProjects\IRonProject\deneme2.dll")

py2phi = phi_c.mc_C
py2phi.restype = None
py2phi.argtypes = [
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
n = 1

y_phi_all = np.empty((3 * n))
z = np.empty((3 * n))
print("z : ", z ,"our : ", y_phi_all)
py2phi(z, y_phi_all)


print("z : ", z ,"our : ", y_phi_all)

