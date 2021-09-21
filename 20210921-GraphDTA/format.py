import pandas as pd
import numpy as np
import json
import pickle

protein = pd.read_csv('./data/balance/protein_sequence.csv', index_col=None, header=None).to_numpy()
smiles = pd.read_csv('./data/balance/smiles.csv', index_col=None, header=None).to_numpy()
mat = pd.read_table('./data/balance/mat_drug_protein.txt', sep=' ', index_col=None, header=None).to_numpy()

protein_dict = {}
need_del = []
for i, item in enumerate(protein):
    [key, value] = item
    if key in protein_dict.keys():
        need_del.append(i)
        continue
    
    protein_dict[key] = value

with open('./data/balance/proteins.txt', 'w') as f:
    json.dump(protein_dict, f)
with open('./data/unbalance/proteins.txt', 'w') as f:
    json.dump(protein_dict, f)

mat = np.delete(mat, need_del, axis=1)
with open('./data/balance/Y', 'wb') as file:
    pickle.dump(mat, file)
with open('./data/unbalance/Y', 'wb') as file:
    pickle.dump(mat, file)

smiles_dict = {}
need_del = []
for i, item in enumerate(smiles):
    [key, value] = item
    if key in smiles_dict.keys():
        need_del.append(i)
        continue
    
    smiles_dict[key] = value

with open('./data/balance/ligands_can.txt', 'w') as f:
    json.dump(smiles_dict, f)
with open('./data/unbalance/ligands_can.txt', 'w') as f:
    json.dump(smiles_dict, f)




