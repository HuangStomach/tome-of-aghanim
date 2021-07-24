import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import Counter

tcga_pred_pd = pd.read_csv('../../Output/3/ccle_a_pred_tcga.csv')
drugs = ['Irinotecan', 'Topptecan'] #tcga_pred_pd.columns[2:]
cancer_types = pd.unique(tcga_pred_pd['Cancer'])
cancer_types = np.delete(cancer_types, [4, 6, 13, 19, 30], axis=0) # UCEC文件错误 暂时删除
sample_type = np.array([x[13:15] for x in tcga_pred_pd.iloc[:, 0].to_numpy()], dtype=str)

fig, axes = plt.subplots(len(drugs), 1, figsize=(12, 12 * len(drugs)))

cancer_drug_ccle = np.array([])
for cancer in cancer_types:
    type_code = "01"
    if cancer == "LAML": type_code = "03"
    if cancer == "SKCM": type_code = "06"

    tmp_pd = tcga_pred_pd.iloc[np.where(sample_type == type_code)]
    tmp_pd = tmp_pd.loc[tmp_pd['Cancer'] == cancer]

    cancer_drug_ccle = np.vstack((cancer_drug_ccle, tmp_pd)) if cancer_drug_ccle.shape[0] > 1 else tmp_pd

for j, drug in enumerate(drugs):
    pos = 0
    for i, cancer in enumerate(cancer_types):
        blca_ccle = cancer_drug_ccle[cancer_drug_ccle[:, 1] == cancer]
        mut_mat = pd.read_table('../../Data/MC3/{}_mc3.txt'.format(cancer))
        sample = np.intersect1d(mut_mat.iloc[:, 0], blca_ccle[:, 0])
        if len(sample) < 50: continue

        mut_samples = mut_mat.loc[mut_mat.iloc[:, 0].isin(sample)]['sample']
        counter = Counter(mut_samples)

        blca_ccle = cancer_drug_ccle[np.where(np.in1d(cancer_drug_ccle[:, 0], sample))[0]]
        middle = np.quantile(blca_ccle[:, j + 2], .75)
        blca_ccle_25 = blca_ccle[blca_ccle[:, j + 2] > middle]
        blca_ccle_75 = blca_ccle[blca_ccle[:, j + 2] <= middle]

        data_25 = np.array([])
        data_75 = np.array([])
        for item in blca_ccle_25:
            data_25 = np.append(data_25, math.log(counter[item[0]]))
        for item in blca_ccle_75:
            data_75 = np.append(data_75, math.log(counter[item[0]]))

        axes[j].boxplot([data_25], positions=[pos + 0.25], widths=.7, patch_artist=True, 
            sym='.', 
            boxprops=dict(facecolor='w', color='lightcoral'), 
            medianprops=dict(color='lightcoral'),
            capprops=dict(color='lightcoral'),
            whiskerprops=dict(color='lightcoral'),
            flierprops=dict(color='lightcoral', markeredgecolor='lightcoral', linewidth='.7'),
        )
        axes[j].boxplot([data_75], positions=[pos + 1.05], widths=.7, patch_artist=True, 
            sym='.',
            boxprops=dict(facecolor='w', color='darkturquoise'),
            medianprops=dict(color='darkturquoise'),
            capprops=dict(color='darkturquoise'),
            whiskerprops=dict(color='darkturquoise'),
            flierprops=dict(color='darkturquoise', markeredgecolor='darkturquoise', linewidth='.7'),
        )
        pos += 2

    axes[j].set_xticks(range(1, len(cancer_types) * 2, 2))
    axes[j].set_xticklabels(cancer_types, rotation = 90)
    axes[j].set_title(drug)
    axes[j].set_ylabel('log(TMB)')
    axes[j].set_facecolor('whitesmoke')

plt.savefig('c.png', bbox_inches='tight')
