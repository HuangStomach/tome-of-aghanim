import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from lib import KNKN
np.set_printoptions(suppress = True)
np.set_printoptions(threshold = np.inf)

class GRMF:
    K = 5
    p = 5
    def __init__(self, max_iter = 1000, eta = 0.003, k = 50, lambda_l = 0.5, lambda_d = 0.001, lambda_t = 0.001):
        '''
        ### Args:
            max_iter: 最大迭代数
            eta: 收敛条件
            k: 奇异值分解个数
            lambda_l: 正则因子
            lambda_d: 正则因子
            lambda_t: 正则因子
        '''
        self.max_iter = max_iter
        self.eta = eta
        self.k = k
        self.lambda_l = lambda_l
        self.lambda_d = lambda_d
        self.lambda_t = lambda_t
    
    def _preprocess(self, K, miu=0.7):
        '''
        ### Args:
            K: 目标的近邻个数
            miu: 衰减因子
        '''
        self.K = K
        Yd = np.zeros((self.n, self.m)) # 存储药物K近邻的加权平均值
        Yt = np.zeros((self.n, self.m)) # 靶点
        Knkn = KNKN(self.K)
        weights = np.zeros(self.K)

        Knkn.fit(self.drug_mat)
        for d in range(self.n):
            (indexes, values) = Knkn.neighbors(d)
            z = np.sum(values)
            for i in range(self.K):
                weights[i] = (miu ** i) * values[i]
                Yd[d] += weights[i] * self.Y[indexes[i]] / z

        Knkn.fit(self.target_mat)
        for t in range(self.m):
            (indexes, values) = Knkn.neighbors(t)
            for i in range(self.K):
                weights[i] = (miu ** i) * values[i]
                Yt[:, t] += weights[i] * self.Y[:, indexes[i]] / z

        for i in range(self.n):
            for j in range(self.m):
                self.Y[i][j] = max(self.Y[i][j], (Yd[i][j] + Yt[i][j]) / 2)

    def _L_hat(self, mat, n) -> np.ndarray:
        '''
        ### Args:
            mat: 矩阵
            n: 元素数目
        '''
        mat_hat = mat.copy().to_numpy()
        D = np.zeros((n, n))
        D_power = np.zeros((n, n))
        Knkn = KNKN(self.p)
        Knkn.fit(mat)
        for i in range(n):
            for j in range(n):
                if i == j:
                    mat_hat[i][j] = 0
                    continue
                
                indexes_i = Knkn.neighbors(i)[0]
                indexes_j = Knkn.neighbors(j)[0]
                n = 0.5
                if j in indexes_i and i in indexes_j: n = 1
                if j not in indexes_i and i not in indexes_j: n = 0
                mat_hat[i][j] = n * mat.iloc[i][j]
            D[i][i] = np.sum(mat_hat[i])
            D_power[i][i] = D[i][i] ** (-1 / 2)

        L = D - mat_hat
        L_hat = D_power @ L @ D_power
        return L_hat
    
    def fit(self, adj_mat, drug_mat, target_mat, preprocess=True, K=5, miu=.7):
        '''
        ### Args:
            adj_mat: 相互作用矩阵
            drug_mat: 药物相似矩阵
            target_mat: 靶点相似矩阵
            preprocess: 是否预处理
            K: 目标的近邻个数
            miu: 衰减因子
        '''
        self.adj_mat = adj_mat
        self.drug_mat = drug_mat
        self.target_mat = target_mat
        (self.n, self.m) = self.adj_mat.shape
        self.Y = self.adj_mat.to_numpy(dtype=float)
        if preprocess: self._preprocess(K, miu)

        [U, S, Vh] = np.linalg.svd(self.Y, full_matrices=False)
        self.k = self.k if min(self.m, self.n) > self.k else min(self.m, self.n)
        U = U[:, :self.k]
        S = np.diag(S[:self.k])
        V = Vh[:self.k, :].transpose()

        s_sqrt = np.sqrt(S)
        A = U @ s_sqrt
        B = V @ s_sqrt

        L_drug_hat = self._L_hat(self.drug_mat, self.n)
        L_target_hat = self._L_hat(self.target_mat, self.m)

        # 迭代计算使其收敛
        for _i in range(self.max_iter):
            A_old = A.copy()
            B_old = B.copy()

            A = (self.Y @ B - self.lambda_d * L_drug_hat @ A) @ np.linalg.inv(B.transpose() @ B + self.lambda_l * np.identity(self.k))
            B = (self.Y.transpose() @ A - self.lambda_t * L_target_hat @ B) @ np.linalg.inv(A.transpose() @ A + self.lambda_l * np.identity(self.k))
            A_diff = A_old - A
            B_diff = B_old - B
            
            # TODO: 求和平方开方做判敛
            diff = max(
                np.concatenate(np.maximum(A_diff, -A_diff)).sum(), 
                np.concatenate(np.maximum(B_diff, -B_diff)).sum()
            )
            if diff < self.eta: break

        self.Y_hat = A @ B.transpose()
        return self

    def predict(self, i, j) -> float:
        '''
        ### Args:
            i: 药物的索引
            j: 靶点的索引
        '''
        return self.Y_hat[i][j]