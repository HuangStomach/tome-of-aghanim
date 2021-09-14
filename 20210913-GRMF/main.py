import numpy as np
import pandas as pd
from lib import KNKN
np.set_printoptions(suppress = True)

mod = 'GPCR'
miu = 0.7
k = 5

adj_mat = pd.read_table('./Data/adj_{}.txt'.format(mod), sep='\t', header=0, index_col=0).T
drug_mat = pd.read_table('./Data/drug_{}.txt'.format(mod), sep='\t', header=0, index_col=0)
target_mat = pd.read_table('./Data/target_{}.txt'.format(mod), sep='\t', header=0, index_col=0)
(n, m) = adj_mat.shape

Yd = np.zeros(n) # 存储药物K近邻的加权平均值
Yt = np.zeros(m) # 靶点
Y = adj_mat.to_numpy(dtype=float)
knkn = KNKN(k)

knkn.fit(drug_mat)
for d in range(n):
    (indexes, values) = knkn.neighbors(d)
    weights = np.zeros(k)
    for i in range(k):
        weights[i] = miu ** (i - 1) * values[i]

    z = np.sum(values)
    adj_neighbors = adj_mat.iloc[indexes].to_numpy()
    Yd[d] = np.dot(weights, adj_neighbors).sum() / z

knkn.fit(target_mat)
for t in range(m):
    (indexes, values) = knkn.neighbors(t)
    weights = np.zeros(k)
    for i in range(k):
        weights[i] = miu ** (i - 1) * values[i]

    z = np.sum(values)
    adj_neighbors = adj_mat.iloc[:, indexes].to_numpy().transpose()
    Yt[t] = np.dot(weights, adj_neighbors).sum() / z

for i in range(n):
    for j in range(m):
        Y[i][j] = max(Y[i][j], (Yd[i] + Yt[j]) / 2)

print(Y)
