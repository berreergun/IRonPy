import sys
import numpy as np
from IRon import phi
import pandas as pd




'''
 Squared Error-Relevance Area (SERA)
 
@description Computes an approximation of the area under the curve described by squared error of predictions for a sequence of subsets with increasing relevance
@param trues Target values from a test set of a given data set. Should be a series and have the same size as the variable preds
@param preds Predicted values given a certain test set of a given data set. Should be a dataframe and have the same size as the variable preds
@param phi_trues Relevance of the values in the parameter trues. Use phi() for more information. Defaults to NULL
@param ph The relevance function providing the data points where the pairs of values-relevance are known. Default is NULL
@param step Relevance intervals between 0 (min) and 1 (max). Default 0.001
@param return_err Boolean to indicate if the errors at each subset of increasing relevance should be returned. Default is FALSE

@export

@return Value for the area under the relevance-squared error curve (SERA)

examples 

dir = os.path.dirname(sys.modules["IRon"].__file__)
path = os.path.join(dir, "data/accel_data.csv")
df = pd.read_csv(path)

# Use one  hot encoder for non-numerical values

att1 = pd.get_dummies(df.Attribute1, prefix='Attribute1')
df = pd.concat([df.drop(['Attribute1'], axis=1), att1], axis=1)
att2 = pd.get_dummies(df.Attribute2, prefix='Attribute2')
df = pd.concat([df.drop(['Attribute2'], axis=1), att2], axis=1)
att5 = pd.get_dummies(df.Attribute5, prefix='Attribute5')
df = pd.concat([df.drop(['Attribute5'], axis=1), att1], axis=1)

#Create target and features
x = df.drop(["acceleration"], axis=1)
y = df["acceleration"]

#Split data as test and train
X_train, X_test, y_train, y_test  = train_test_split(x,y, test_size=0.3, random_state=7)

#Calculate phi
phit = phi.phi_control(y_train)
phis = phi.phi(y_test, phi_parms=phit)

#ML Part

#Random forest part
from sklearn.ensemble import RandomForestRegressor
regressor = RandomForestRegressor(n_estimators = 40, random_state = 0)
regressor.fit(X_train, y_train)
y_pred_random_forest = regressor.predict(X_test)

# SVM Regression Part
from sklearn.svm import SVR
svr = SVR()
svr.fit(X_train, y_train)
y_pred_SVR = svr.predict(X_test)


y_pred = pd.DataFrame(
        {'RFPredictions': y_pred_random_forest,
         'SVRPredictions': y_pred_SVR,
         })


#Calculate SERA
from IRon import sera
SERA = sera.sera(preds=y_pred,phi_trues=phis,trues=y_test)


plt.bar( SERA.keys(),SERA, align='center', alpha=0.5)
plt.title('SERA')
plt.show()

#alternative
SERA = sera.sera(preds=y_pred,ph=phit,trues=y_test)

'''



def sera (trues, preds, phi_trues = None, ph = None, step = 0.001,return_err = False) :
    if not isinstance(preds, pd.DataFrame):
        preds = pd.DataFrame(preds)

    if ((phi_trues is None) & (ph is None)):
        sys.exit("You need to input either the parameter phi_trues or ph.")

    if phi_trues is None :
       phi_trues = phi.phi(trues, ph)

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

    if return_err :
       return {"sera":res, "errors":errors, "thrs" :th}
    else:
       return res







