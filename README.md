# Delta_NMR

Delta NMR is a KRR ML model for prediciton of chemical tensors at the wb97xd | 6-311g(d,p) level of theory. Trained on Dataset 4, Delta NMR takes as an input log files with results of chemical tensors at cheap wb97xd | 3-21G level of theory. This project was developed by Krystof Chrappova as part of Msci final year project.
![image](https://user-images.githubusercontent.com/76857765/144690608-f0ef11e6-ac4e-4d4a-85a0-7d1a5beeef8a.png)


Known limitaions:
1. Molecular size - any molecule that has more than 80 features created with ACSF representation will raise an error
2. Conjugation - cutoff function is set to 3 Angstrom. This is molecular representation limitation and shifts for conjugated systems will not be accurate.

## Dependencies
* scikit-learn              0.24.2
* scikit-optimize           0.9.0
* dscribe                   1.1.0
* ase                       3.22.0
* pandas                    1.3.2
* numpy                     1.21.2

## Usage
Delta_NMR works as follows:
1. run Gaussian calculation with wb97xd | 3-21G level of theory and create `.sdf` files from the `.log` files using mol_translator (without scaling the tensors)
3. create atoms dataframe with calculation results by running make_df.py (adjust path to mol_translator and your sdf caluclations results)
5. using make_delta_inp.py prepare the input dataframe for the delta ML model.
6. using delta_predict.py get wb97xd | 6-311g(d,p) level of theory chemical tensors

## Carbon and Proton Tensor Prediction Error (test errors)

| Atom |  RMSE  |  MAE  | MAXE |
|------|--------|-------|------|
| H    |  0.155 | 0.111 | 1.28 |
| C    |  1.48  | 1.09  | 10.8 |



## Retraining
Models in `trained_models` were trained on dataset 4 and testes on dataset 3  as shown in the flowchart. If you wish to retrain the model you can do so using tain_test_delta.py. Here you need to change the path to high level and low level claulcation reults (dataframes created with `make_delta_inp.py`)
