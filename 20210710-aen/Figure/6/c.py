import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import Counter


tcga_pred_pd = pd.read_csv('../../Output/3/ccle_a_pred_tcga.csv')
drugs = tcga_pred_pd.columns[2:]
cancer_types = pd.unique(tcga_pred_pd['Cancer'])
sample_type = np.array([x[13:15] for x in tcga_pred_pd.iloc[:, 0].to_numpy()], dtype=str)

fig, axes = plt.subplots(1, len(drugs))

cancer_drug_ccle = np.array([])
for cancer in cancer_types:
    if cancer == 'UCEC': continue

    type_code = "01"
    if cancer == "LAML": type_code = "03"
    if cancer == "SKCM": type_code = "06"

    tmp_pd = tcga_pred_pd.iloc[np.where(sample_type == type_code)]
    tmp_pd = tmp_pd.loc[tmp_pd['Cancer'] == cancer]

    cancer_drug_ccle = np.vstack((cancer_drug_ccle, tmp_pd)) if cancer_drug_ccle.shape[0] > 1 else tmp_pd

for i, drug in enumerate(drugs):
    data = np.array([])
    for cancer in cancer_types:
        if cancer == 'UCEC': continue

        blca_ccle = cancer_drug_ccle[cancer_drug_ccle[:, 1] == cancer]
        mut_mat = pd.read_table('../../Data/MC3/{}_mc3.txt'.format(cancer))
        sample = np.intersect1d(mut_mat.iloc[:, 0], blca_ccle[:, 0])

        if len(sample) < 50: continue

        mut_samples = mut_mat.loc[mut_mat.iloc[:, 0].isin(sample)]['sample']
        counter = Counter(mut_samples)

        blca_ccle = cancer_drug_ccle[np.where(np.in1d(cancer_drug_ccle[:, 0], sample))[0]]