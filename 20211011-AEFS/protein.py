import numpy as np
import pandas as pd
import torch
import torch.nn.utils.rnn as rnn_utils
np.set_printoptions(suppress=True, linewidth=np.nan)

def collate_fn(data):
    data.sort(key=lambda x: len(x), reverse=True)
    data = rnn_utils.pad_sequence(data, batch_first=True, padding_value=0)
    return data

def transform():
    train_prots = np.loadtxt("./datasets/DTINet/protein_sequence.csv", dtype=str, delimiter=',')[:, 1]
    blosum62 = pd.read_table("./datasets/DTINet/blosum62.txt", header=0, index_col=0, sep=' ', dtype=str)
    proteins = []
    for protein in train_prots:
        mat = np.zeros((len(protein), blosum62.shape[0]))
        for i, char in enumerate(protein):
            mat[i] = blosum62.loc[char].to_numpy()
        proteins.append(mat)
    