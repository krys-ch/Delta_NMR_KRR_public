import pandas as pd
import numpy as np
import math
from skopt import BayesSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn import metrics
from sklearn.metrics import max_error, mean_absolute_error, mean_squared_error, median_absolute_error
from joblib import dump, load
import pickle

#high level and low level datasets



dst_3_hl = pd.read_csv(r'..........')
dst_3_ll = pd.read_csv(r'..........')
dst_4_hl = pd.read_csv(r'..........')
dst_4_ll = pd.read_csv(r'..........')



zeros = []
for element in np.arange((len(dst_3_ll.columns)-12), 80, 1).tolist():
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



krr = KernelRidge(alpha = 5.2082289494033395e-08, degree = 8, gamma= 9.908363944271308e-05, kernel = 'poly')
krr.fit(X_train, y_train)

pickle.dump(krr,open("IMPACT_H.pkl","wb"))

#scaling H -1.0719, 32.1254
#scaling C -1.0399, 187.136
pred_train =[]
pred_test = []
for prediction in krr.predict(X_train):
    pred_train.append(float(prediction))
for prediction in krr.predict(X_test):
    pred_test.append(float(prediction))
print('y_test', y_test)
print('pred_test', pred_test)
true_test = (y_test - 32.1254)/-1.0719
predicted_test = pd.DataFrame((np.array(pred_test)-32.1254)/-1.0719)
true_train = (y_train - 32.1254)/-1.0719
predicted_train = pd.DataFrame((np.array(pred_train)-32.1254)/-1.0719)



print('getting stats for train and test....')
print('len true_test', len(true_test))
print('len predicted_test', len(predicted_test))

a = dst_3_ll[['molecule_name', 'atom_index', 'typestr', 'typeint', 'x', 'y', 'z','conn']]
a['IMPACT_shift_ppm'] = predicted_test
a['6_311G_shift_ppm'] = true_test
a.to_csv('IMPACT_output_test_H.csv')

b = dst_4_ll[['molecule_name', 'atom_index', 'typestr', 'typeint', 'x', 'y', 'z','conn']]
b['IMPACT_shift_ppm'] = predicted_train
b['6_311G_shift_ppm'] = true_train
b.to_csv('IMPACT_output_train_H.csv')



true_test.to_csv('labels_test_H.csv')
predicted_test.to_csv('IMPACT_test_H.csv')
true_train.to_csv('labels_train_H.csv')
predicted_train.to_csv('IMPACT_train_H.csv')



r2_test = metrics.r2_score(true_test, predicted_test)
r2_train = metrics.r2_score(true_train, predicted_train)
MAE_train = mean_absolute_error(true_train, predicted_train)

MAE_test = mean_absolute_error(true_test,predicted_test)

RMSE_train = (mean_squared_error(true_train, predicted_train, squared=False))

RMSE_test = (mean_squared_error(true_test, predicted_test, squared=False))

MAXE_train = max_error(true_train, predicted_train)

MAXE_test = max_error(true_test, predicted_test)



median_error_train = median_absolute_error(true_train, predicted_train)
median_error_test = median_absolute_error(true_test, predicted_test)
#stats for the model
print('r2_test: ', r2_test)
print('r2_train: ', r2_train)
print('MAE_train: ', MAE_train)
print('MAE_test: ', MAE_test)
print('RMSE_train: ', RMSE_train)
print('RMSE_test: ', RMSE_test)
print('MAXE_train: ', MAXE_train)
print('MAXE_test: ', MAXE_test)
print('median_train_abs_error', median_error_train)
print('median_test_abs_error', median_error_test)
