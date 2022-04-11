import pandas as pd
import numpy as np
import math
from skopt import BayesSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn import metrics
from sklearn.metrics import max_error, mean_absolute_error, mean_squared_error

#high level and low level datasets
dst_3_hl = pd.read_csv(r'C:\model_input_data\dst_3_rcut_all\dst_3_hl_H_3.csv')
dst_3_ll = pd.read_csv(r'C:\model_input_data\dst_3_rcut_all\dst_3_ll_H_3.csv')
dst_4_hl = pd.read_csv(r'C:\model_input_data\dst_4_rcut_all\dst_4_hl_H_3.csv')
dst_4_ll = pd.read_csv(r'C:\model_input_data\dst_4_rcut_all\dst_4_ll_H_3.csv')

#adding zeros to match number of features of dataset on which model will be trained
zeros = []
for element in np.arange(56, 80, 1).tolist():
    a = str(element)
    zeros.append(a)
dst_3_ll[zeros] = 0
dst_3_hl[zeros] = 0


features = np.arange(0, 80, 1).tolist()
features.append('shift')
f = []
for x in features:
    a = str(x)
    f.append(a)
labels = ['shift']

for feature in dst_4_ll.columns:
    str(feature)

X_train = dst_4_ll[f]
X_test = dst_3_ll[f]
y_train = dst_4_hl[labels]
y_test = dst_3_hl[labels]


print('optimizing hyperparameters')
# log-uniform: understand as search over p = exp(x) by varying x
opt = BayesSearchCV(
        KernelRidge(),
    {  'alpha': (1e-6, 1e+1, 'log-uniform'),
        'gamma': (1e-6, 1e+1, 'log-uniform'),
        'degree': (1, 8),  # integer valued parameter
        'kernel': ['linear', 'poly', 'rbf', 'sigmoid'],  # categorical parameter
     },
    n_iter=32,
    cv=3
)

opt.fit(X_train, y_train)

print("best params: %s" % str(opt.best_params_))


krr = KernelRidge(opt.best_params_)
krr.fit(X_train, y_train)

dump(krr, open('delta_model_H.pkl', 'wb'))

print('getting stats for train and test....')

r2_test = metrics.r2_score(y_test,
                                krr.predict(X_test))
r2_train = metrics.r2_score(y_train,
                                 krr.predict(X_train))
MAE_train = (mean_absolute_error(y_train,
                                     krr.predict(X_train)))

MAE_test = (mean_absolute_error(y_test,
                                    krr.predict(X_test)))

RMSE_train = (mean_squared_error(y_train,
                                     krr.predict(X_train), squared=False))

RMSE_test = (mean_squared_error(y_test,
                                    krr.predict(X_test), squared=False))

MAXE_train = (max_error(y_train,
                            krr.predict(X_train)))

MAXE_test = (max_error(y_test,
                           krr.predict(X_test)))

#stats for the model
print('r2_test: ', r2_test)
print('r2_train: ', r2_train)
print('MAE_train: ', MAE_train)
print('MAE_test: ', MAE_test)
print('RMSE_train: ', RMSE_train)
print('RMSE_test: ', RMSE_test)
print('MAXE_train: ', MAXE_train)
print('MAXE_test: ', MAXE_test)

