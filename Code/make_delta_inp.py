'''
Checklist for running this:
- set path for the sdf files
- set path for correct dataframe (.pkl)
- change name of outfile
'''
from dscribe.descriptors import ACSF
from tqdm import tqdm
import pandas as pd
import os
from ase.io import read
import pickle

path_sdf = r"C:\\Strychnine"
atoms_dfs = os.listdir(r"C:\\Strychnine")

acsfs = []
mol_name = []
bad_files = []
for root, directories, files in os.walk(path_sdf, topdown=False):
    print('Creating ACSF representation for each sdf file...')
    for name in tqdm(files):
        if name.endswith('.sdf'):
            try:
                print(name)
                mol_path = os.path.join(root, name)
                structure = read(str(mol_path))
                species = set()
                species.update(structure.get_chemical_symbols())
                # cmat = CoulombMatrix(n_atoms_max=len(structure))
                acsf = ACSF(
                    species=species,
                    rcut=6.0,
                    g2_params=[[1, 1], [1, 2], [1, 3]],
                    g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, 1]],
                )
                acsf_testmol = acsf.create(structure)
                acsfs.append(acsf_testmol)
                mol_name.append(name)
            except ValueError:
                print('bad file')
                print(name)
                bad_files.append(name)

all_acsfs = pd.DataFrame((zip(mol_name, acsfs)), columns=['mol_name', 'ACSF'])


dataframes = []
for pic in atoms_dfs:
    if pic.endswith('.pkl'):
        print(pic)
        pickle_off = open(r"C:\\Strychnine\\" + str(pic), 'rb')
        df = pickle.load(pickle_off)
        dataframes.append(df)
    else:
        continue

print(dataframes)

# sort df by moelcule name and atom index to ensure correct ACSF is appended
df_atoms = dataframes[0].sort_values(['molecule_name', 'atom_index'])
df_atoms.reset_index()
all_acsfs.sort_values(['mol_name'], inplace=True)


# adding ACSF for each atom
ACSF = []
for lists in all_acsfs.ACSF:
    for value in lists:
        ACSF.append(value)


df_atoms['ACSF'] = ACSF
df_atoms.reset_index()

# expanding each ACSF array into separate columns, adding 0 for each NaN in the dataframe
split_df = pd.DataFrame(df_atoms['ACSF'].tolist())
last_df = split_df.join(df_atoms.reset_index())

clean = last_df.fillna(0)
clean.reset_index()
# creating dataframe for each atom type
clean_H = clean.loc[clean['typeint'] == 1]
clean_C = clean.loc[clean['typeint'] == 6]


clean_H.to_csv('delta_input_H.csv')
clean_C.to_csv('delta_input_C.csv')



if len(bad_files) == 0:
    print('There are no bad files, wohooo!')
else:
    print('There are ', len(bad_files), ' bad files!! This is because some sdf files have no space between the number of atoms and number of bonds - this only happens with molecules with more than 100 bonds. You have to fix the sdf files.')

a = clean['typeint'].unique()
print("Atoms in calculation: ", sorted(a))
