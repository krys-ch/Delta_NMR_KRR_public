import glob
from tqdm import tqdm
import sys
from mol_translator.aemol import aemol
from mol_translator.imp_converter import dataframe_prep as df_prep
from mol_translator.imp_converter import dataframe_write as df_write

#path to mol translator
sys.path.append('/mnt/storage/home/ep18192/opt/bg_code/mol_translator-master/')
files = glob.glob('/mnt/storage/home/ep18192/project_stuff/autoenrich-master/OPT/*.sdf')

amols = []
for file in tqdm(files):
    p = file.split('OPT/')[-1].split('.')[0]
    print(p)
    amol = aemol(p)
    amol.from_file(file, ftype='sdf')
    amol.prop_fromfile(file, 'g09', 'scf')
    mol = df_prep.prep_mol_nmr(amol, nmr_file = str(file), nmr_type = 'nmredata')
    amols.append(mol)

print(len(amols))

df_write.make_atom_df(amols, progress=True, write=True)