import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.utils.rnn as rnn_utils
from sklearn.metrics import jaccard_score
np.set_printoptions(suppress=True, linewidth=np.nan)

hidden_size = 128
num_layers = 2

def proteins():
    train_prots = np.loadtxt("./datasets/DTINet/protein_sequence.csv", dtype=str, delimiter=',')[:, 1]
    blosum62 = pd.read_table("./datasets/DTINet/blosum62.txt", header=0, index_col=0, sep=' ', dtype=str)

    proteins = []
    seq_lens = []
    for protein in train_prots:
        mat = []
        for char in protein:
            if char not in blosum62.columns: continue
            mat.append(blosum62.loc[char].to_numpy(np.int32))
        seq_lens.append(len(mat))
        mat = np.array(mat)
        proteins.append(torch.from_numpy(mat))
    proteins = torch.nn.utils.rnn.pad_sequence(proteins, batch_first=True).to(torch.float32)
    seq_lens = torch.tensor(seq_lens)
    
    output = None
    lstm = nn.LSTM(blosum62.shape[0], hidden_size, 2, bidirectional=True, batch_first=True)
    sorted_lens, indices = seq_lens.sort(descending=True)
    _, un_idx = torch.sort(indices, dim=0)
    sorted_proteins = torch.from_numpy(proteins.numpy()[indices])
    proteins_pack = rnn_utils.pack_padded_sequence(sorted_proteins, sorted_lens, batch_first=True)


    output, (h, c) = lstm(proteins_pack)
    out, _ = rnn_utils.pad_packed_sequence(output, batch_first=True)
    out = torch.index_select(out, 0, un_idx)

    embeds = []
    for o, length in zip(out, seq_lens):
        embed = o[length - 1, :hidden_size] + o[length - 1, hidden_size:]
        embeds.append(embed.detach().numpy())
    embeds = np.array(embeds)
    np.savetxt('./datasets/DTINet/protein_embeds.csv', embeds, delimiter=',')

def diseases(dataset):
    pdi = dataset.pdi()
    disease_num = pdi.shape[1]

    disease_A = np.zeros((disease_num, disease_num), dtype=float)
    for i in range(disease_num):
        for j in range(disease_num):
            disease_A[i, j] = jaccard_score(pdi[:, i], pdi[:, j])
    np.savetxt('./datasets/DTINet/Similarity_Matrix_Diseases.txt.txt', disease_A, delimiter=' ')

if __name__=='__main__':
    proteins()
    diseases()