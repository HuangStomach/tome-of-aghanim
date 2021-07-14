import os
import joblib
import pandas as pd

# 保存tcga样本到癌症类型的映射
dirs = os.listdir('./Data/TCGA/')
tcga_to_cancer = dict()

for cancer in dirs:
    if cancer == '.gitkeep': continue
    tcga = pd.read_table("./Data/TCGA/{}/HiSeqV2".format(cancer))
    tcga_rpkm = tcga.columns.to_numpy(dtype=str)
    for rpkm in tcga_rpkm[1:]:
        tcga_to_cancer[rpkm.replace('-', '.')] = cancer

joblib.dump(tcga_to_cancer, './Output/1/tcga_to_cancer.joblib')  

# 保存ccle样本到细胞系组织的映射
ccle_latent = pd.read_table('./Output/1/1.CCLE_latent.tsv', index_col = 0)
ccle_to_tissue = dict()

for rowname in ccle_latent.index.values:
    tissue = rowname.split('.')[0]
    tissue = tissue.split('_')[1:]
    ccle_to_tissue[rowname] = '_'.join(tissue)

# print(unique(list(ccle_to_tissue.values())))
joblib.dump(ccle_to_tissue, './Output/1/ccle_to_tissue.joblib')  
