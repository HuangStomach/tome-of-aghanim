import pandas as pd
import numpy as np

ddict = {}
drug_dict = np.loadtxt('./data/DTINet/drug_dict_map.txt', dtype=str, delimiter=':')
for row in drug_dict:
    ddict[row[1].upper()] = row[0]

drugs_dict = {}
drugs = np.loadtxt('./data/DTINet/drug.txt', dtype=str, delimiter='\n')
for i, drug in enumerate(drugs):
    drugs_dict[drug] = i

proteins_dict = {}
proteins = np.loadtxt('./data/DTINet/protein.txt', dtype=str, delimiter='\n')
for i, protein in enumerate(proteins):
    proteins_dict[protein] = i

dtinet_dti = np.loadtxt('./data/DTINet/mat_drug_protein.txt', dtype=int, delimiter=' ')
dtinet_dti_c = dtinet_dti.copy()

dc_dti = pd.read_csv('./data/drugcentral_DTI.csv')
for index, row in dc_dti.iterrows():
    upper = row['DRUG_NAME'].upper()
    if upper not in ddict.keys(): continue
    drug = ddict[upper]
    if drug not in drugs_dict.keys(): continue
    d_index = drugs_dict[drug]
    proteins = row['ACCESSION'].split('|')
    for protein in proteins:
        if protein not in proteins_dict.keys(): continue
        p_index = proteins_dict[protein]
        dtinet_dti_c[d_index][p_index] = 1

np.savetxt('./data/DTINet/mat_drug_protein_s.txt', dtinet_dti_c, fmt='%d', delimiter=' ')
