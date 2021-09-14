import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from lib import KNKN
np.set_printoptions(suppress = True)

mod = 'GPCR'
miu = 0.7
K = 5

adj_mat = pd.read_table('./Data/adj_{}.txt'.format(mod), sep='\t', header=0, index_col=0).T
drug_mat = pd.read_table('./Data/drug_{}.txt'.format(mod), sep='\t', header=0, index_col=0)
target_mat = pd.read_table('./Data/target_{}.txt'.format(mod), sep='\t', header=0, index_col=0)
(n, m) = adj_mat.shape

## WKNKN
Yd = np.zeros((n, m)) # 存储药物K近邻的加权平均值
Yt = np.zeros((n, m)) # 靶点
Y = adj_mat.to_numpy(dtype=float)
KK = KNKN(K)
weights = np.zeros(K)

KK.fit(drug_mat)
for d in range(n):
    (indexes, values) = KK.neighbors(d)
    z = np.sum(values)
    for i in range(K):
        weights[i] = (miu ** i) * values[i]
        Yd[d] += weights[i] * Y[indexes[i]] / z

KK.fit(target_mat)
for t in range(m):
    (indexes, values) = KK.neighbors(t)
    for i in range(K):
        weights[i] = (miu ** i) * values[i]
        Yt[:, t] += weights[i] * Y[:, indexes[i]] / z

for i in range(n):
    for j in range(m):
        Y[i][j] = max(Y[i][j], (Yd[i][j] + Yt[i][j]) / 2)

# GRMF
k = 50
lambda_l = 0.5
lambda_d = 0.001
lambda_t = 0.001
svd = TruncatedSVD(n_components=k, random_state=42)
svd.fit(Y)
Y = svd.transform(Y)
[u, s, vh] = np.linalg.svd(Y, full_matrices=False)
print(s)