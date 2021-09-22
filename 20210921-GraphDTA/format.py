import pandas as pd
import numpy as np
import json
import pickle

protein = pd.read_csv('./source/protein_sequence.csv', index_col=None, header=None).to_numpy()
smiles = pd.read_csv('./source/smiles.csv', index_col=None, header=None).to_numpy()
mat = pd.read_table('./source/mat_drug_protein.txt', sep=' ', index_col=None, header=None).to_numpy()

protein_dict = {}
need_del_col = []
for i, item in enumerate(protein):
    [key, value] = item
    
    if key in protein_dict.keys():
        need_del_col.append(i)
        continue
    
    protein_dict[key] = value

with open('./data/balance/proteins.txt', 'w') as f:
    json.dump(protein_dict, f)
with open('./data/unbalance/proteins.txt', 'w') as f:
    json.dump(protein_dict, f)

mat = np.delete(mat, need_del_col, axis=1)

ignore = ['DB01356', 'DB01378']
smiles_dict = {}
need_del_row = []
for i, item in enumerate(smiles):
    [key, value] = item
    if key in ignore or key in smiles_dict.keys():
        need_del_row.append(i)
        continue
    
    smiles_dict[key] = value
mat = np.delete(mat, need_del_row, axis=0)
with open('./data/balance/Y', 'wb') as file:
    pickle.dump(mat, file)
with open('./data/unbalance/Y', 'wb') as file:
    pickle.dump(mat, file)

with open('./data/balance/ligands_can.txt', 'w') as f:
    json.dump(smiles_dict, f)
with open('./data/unbalance/ligands_can.txt', 'w') as f:
    json.dump(smiles_dict, f)

# 平衡训练集
ones = []
zeros = []
for i, v in enumerate(mat.flatten()):
    if v == 0: zeros.append(i)
    else: ones.append(i)
ones = np.array(ones)
zeros = np.array(zeros)

np.random.shuffle(ones)
np.random.shuffle(zeros)

splits = np.array(np.array_split(ones, 5))
train = splits[0:4]
test = np.array(splits[4])

splits = np.array_split(zeros[:len(ones)], 5)
train[0] = np.hstack((train[0], splits[0])).tolist()
train[1] = np.hstack((train[1], splits[1])).tolist()
train[2] = np.hstack((train[2], splits[2])).tolist()
train[3] = np.hstack((train[3], splits[3])).tolist()
test = np.append(test, zeros[int(.2 * len(zeros)):int(.4 * len(zeros))]).tolist()

#train = train.tolist()
with open('./data/balance/folds/train_fold_setting1.txt', 'w') as f:
    json.dump(train.tolist(), f)
with open('./data/balance/folds/test_fold_setting1.txt', 'w') as f:
    json.dump(test, f)

# 不平衡训练集
ones = []
zeros = []
for i, v in enumerate(mat.flatten()):
    if v == 0: zeros.append(i)
    else: ones.append(i)
ones = np.array(ones)
zeros = np.array(zeros)

np.random.shuffle(ones)
np.random.shuffle(zeros)

splits = np.array(np.array_split(ones, 5))
train = splits[0:4]
test = np.array(splits[4])

splits = np.array_split(zeros, 5)
train[0] = np.hstack((train[0], splits[0])).tolist()
train[1] = np.hstack((train[1], splits[1])).tolist()
train[2] = np.hstack((train[2], splits[2])).tolist()
train[3] = np.hstack((train[3], splits[3])).tolist()
test = np.append(test, splits[4]).tolist()

with open('./data/unbalance/folds/train_fold_setting1.txt', 'w') as f:
    json.dump(train.tolist(), f)
with open('./data/unbalance/folds/test_fold_setting1.txt', 'w') as f:
    json.dump(test, f)
