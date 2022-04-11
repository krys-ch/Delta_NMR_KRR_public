import pandas as pd
import numpy as np
from pickle import load
from sklearn.kernel_ridge import  KernelRidge
tensor_ll = pd.read_csv('delta_input_C.csv')
zeros = []

for element in np.arange((len(tensor_ll.columns)-12), 80, 1).tolist():
    a = str(element)
    zeros.append(a)
tensor_ll[zeros] = 0

features = np.arange(0, 80, 1).tolist()
features.append('shift')
f = []
for x in features:
    a = str(x)
    f.append(a)
labels = ['shift']

model = load(open('trained_models/delta_model_C.pkl', 'rb'))

tensor_hl = model.predict(tensor_ll[f]).tolist()
tensor_ll['predicted_tensor'] = tensor_hl
result_df = tensor_ll[['molecule_name', 'atom_index', 'typeint', 'predicted_tensor', 'x', 'y', 'z']]
result_df.to_csv('predicted_tensors_C.csv')
