import os
import numpy as np
import pandas as pd

dirs = os.listdir('./Data/TCGA/')
tcga_ss_mat = dict()

for cancer in dirs:
    if cancer == '.gitkeep': continue
    print(cancer)
    tcga = pd.read_table("./Data/TCGA/{}/HiSeqV2".format(cancer))
    tcga_rpkm = tcga.columns.to_numpy(dtype=str)
    for rpkm in tcga_rpkm[1:]:
        tcga_ss_mat[rpkm.replace('-', '.')] = cancer

np.save('./Output/1/tcga_ss_mat.npy', tcga_ss_mat)