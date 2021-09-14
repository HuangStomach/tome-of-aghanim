import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from lib import KNKN
np.set_printoptions(suppress = True)
np.set_printoptions(threshold = np.inf)

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
Knkn = KNKN(K)
weights = np.zeros(K)

Knkn.fit(drug_mat)
for d in range(n):
    (indexes, values) = Knkn.neighbors(d)
    z = np.sum(values)
    for i in range(K):
        weights[i] = (miu ** i) * values[i]
        Yd[d] += weights[i] * Y[indexes[i]] / z

Knkn.fit(target_mat)
for t in range(m):
    (indexes, values) = Knkn.neighbors(t)
    for i in range(K):
        weights[i] = (miu ** i) * values[i]
        Yt[:, t] += weights[i] * Y[:, indexes[i]] / z

for i in range(n):
    for j in range(m):
        Y[i][j] = max(Y[i][j], (Yd[i][j] + Yt[i][j]) / 2)

# GRMF
k = 50
# 先做奇异值分解，取前k个奇异值对应的特征向量
[U, S, Vh] = np.linalg.svd(Y, full_matrices=False)
U = U[:, :k]
S = np.diag(S[:k])
V = Vh[:k, :].transpose()

s_sqrt = np.sqrt(S)
A = U @ s_sqrt
B = V @ s_sqrt

# 计算拉格朗日乘子法所需参数
lambda_l = 0.5
lambda_d = 0.001
lambda_t = 0.001
drug_mat_hat = drug_mat.copy().to_numpy()
D_d = np.zeros((n, n))
D_d_power = np.zeros((n, n))
Knkn.fit(drug_mat)
for i in range(n):
    for j in range(n):
        if i == j:
            drug_mat_hat[i][j] = 0
            continue
        
        indexes_i = Knkn.neighbors(i)[0]
        indexes_j = Knkn.neighbors(j)[0]
        n = 0.5
        if j in indexes_i and i in indexes_j: n = 1
        if j not in indexes_i and i not in indexes_j: n = 0
        drug_mat_hat[i][j] = n * drug_mat.iloc[i][j]
    D_d[i][i] = np.sum(drug_mat_hat[i])
    D_d_power[i][i] = D_d[i][i] ** (-1 / 2)

target_mat_hat = target_mat.copy().to_numpy()
D_t = np.zeros((m, m))
D_t_power = np.zeros((m, m))
Knkn.fit(target_mat)
for i in range(m):
    for j in range(m):
        if i == j:
            target_mat_hat[i][j] = 0
            continue
        
        indexes_i = Knkn.neighbors(i)[0]
        indexes_j = Knkn.neighbors(j)[0]
        n = 0.5
        if j in indexes_i and i in indexes_j: n = 1
        if j not in indexes_i and i not in indexes_j: n = 0
        target_mat_hat[i][j] = n * target_mat.iloc[i][j]
    D_t[i][i] = np.sum(target_mat_hat[i])
    D_t_power[i][i] = D_t[i][i] ** (-1 / 2)

L_d = D_d - drug_mat_hat
L_t = D_t - target_mat_hat
Ld_hat = D_d_power @ L_d @ D_d_power
Lt_hat = D_t_power @ L_t @ D_t_power

# 迭代计算使其收敛
max_iter = 1000
for _i in range(max_iter):
    A_old = A.copy()
    B_old = B.copy()

    A = (Y @ B - lambda_d * Ld_hat @ A) @ np.linalg.inv(B.transpose() @ B + lambda_l * np.identity(k))
    B = (Y.transpose() @ A - lambda_t * Lt_hat @ B) @ np.linalg.inv(A.transpose() @ A + lambda_l * np.identity(k))
    A_diff = A_old - A
    B_diff = B_old - B
    
    diff = max(
        np.concatenate(np.maximum(A_diff, -A_diff)).sum(), 
        np.concatenate(np.maximum(B_diff, -B_diff)).sum()
    )
    if diff < 0.003: break

Y_hat = A @ B.transpose()
print(Y_hat.shape)
