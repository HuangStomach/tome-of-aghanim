import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.utils.rnn as rnn_utils
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import jaccard_score
# from transformers import AutoTokenizer, AutoModelForMaskedLM

from urllib import request
from time import sleep

np.set_printoptions(suppress=True, linewidth=np.nan)

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
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
    
    lstm = nn.LSTM(blosum62.shape[0], hidden_size, 2, bidirectional=True, batch_first=True).to(device)
    sorted_lens, indices = seq_lens.sort(descending=True)
    _, un_idx = torch.sort(indices, dim=0)
    sorted_proteins = torch.from_numpy(proteins.numpy()[indices])
    proteins_pack = rnn_utils.pack_padded_sequence(sorted_proteins, sorted_lens, batch_first=True).to(device)

    output, (h, c) = lstm(proteins_pack)
    out, _ = rnn_utils.pad_packed_sequence(output, batch_first=True)
    out = torch.index_select(out, 0, un_idx)

    embeds = []
    for o, length in zip(out, seq_lens):
        embed = o[length - 1, :hidden_size] + o[length - 1, hidden_size:]
        embeds.append(embed.detach().numpy())
    embeds = np.array(embeds)
    np.savetxt('./datasets/DTINet/protein_embeds.csv', embeds, delimiter=',')

def smiles():
    seqs = []
    drug_dict = np.loadtxt('./data/drug.txt', dtype=str, delimiter='\n')
    drug_url = 'https://go.drugbank.com/structures/small_molecule_drugs/{}.smiles'

    for drug in drug_dict:
        sleep(1)
        try:
            req = request.Request(drug_url.format(drug), headers={
                'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.110 Safari/537.36'
            })
            data = request.urlopen(req).read()
            seqs.append([drug, data])
            print(drug, 'OK')

        except Exception as e:
            print(drug, e)
            seqs.append([drug, 'ERROR'])

    np.savetxt('./data/drug_smiles.csv', seqs, fmt='%s', delimiter=',')

def ecfps():
    seqs = []
    radius = 6
    length = 4096

    drugs = np.loadtxt('./data/DTINet/drug_smiles.csv', delimiter=',', dtype=str, comments=None)
    for drug in drugs:
        try:
            name, smiles = drug
            mol = Chem.MolFromSmiles(smiles)
            seqs.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length).ToList())
        except Exception as e:
            print(drug, e)

    np.savetxt('./data/DTINet/drug_ecfps{}.txt'.format(radius * 2), seqs, fmt='%s', delimiter=',')

def drug_token():
    tokenizer = AutoTokenizer.from_pretrained("seyonec/SMILES_tokenized_PubChem_shard00_160k")


def drug_sim():
    # drugs_pubchem = pd.read_table('./data/LRSSL/drug_pubchem_mat.txt', sep='\t', index_col=0)
    # drugs_domain = pd.read_table('./data/LRSSL/drug_target_domain_mat.txt', sep='\t', index_col=0)
    # drugs_go = pd.read_table('./data/LRSSL/drug_target_go_mat.txt', sep='\t', index_col=0)
    # length = drugs_pubchem.shape[0]
    # sim = np.zeros((length, length))

    # for i in range(length):
    #     for j in range(length):
    #         drugA = np.concatenate((drugs_pubchem.iloc[i], drugs_domain.iloc[i], drugs_go.iloc[i]), axis=1)
    #         drugB = np.concatenate((drugs_pubchem.iloc[j], drugs_domain.iloc[j], drugs_go.iloc[j]), axis=1)
    #         sim[i][j] = jaccard_score(drugA.to_numpy(), drugB.to_numpy())

    # np.savetxt('./data/LRSSL/drug_sim.txt', sim, fmt='%s', delimiter=' ')

    # drugs_protein = np.loadtxt('./data/DTINet/mat_drug_protein_s.txt', dtype=int, delimiter=' ')
    # drugs_disease = np.loadtxt('./data/DTINet/mat_drug_disease.txt', dtype=int, delimiter=' ')
    drugs_ecfps = np.loadtxt('./data/DTINet/drug_ecfps12.txt', dtype=int, delimiter=',')
    length = drugs_ecfps.shape[0]
    sim = np.zeros((length, length))

    for i in range(length):
        for j in range(length):
            drugA = drugs_ecfps[i]
            drugB = drugs_ecfps[j]
            sim[i][j] = jaccard_score(drugA, drugB)

    np.savetxt('./data/DTINet/Similarity_Matrix_Drugs_s.txt', sim, fmt='%s', delimiter='    ')

if __name__=='__main__':
    # proteins()
    ecfps()
    # drug_sim()
