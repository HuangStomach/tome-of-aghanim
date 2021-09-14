import math
import pandas as pd
import numpy as np
from scipy.io import loadmat 
import lib
np.set_printoptions(suppress = True)

drugs_sim = pd.read_table('./Data/DrugSimMat', sep=' ', header=None).to_numpy() # 药物相似程度矩阵
diseases_sim = pd.read_table('./Data/DiseaseSimMat', sep=' ', header=None).to_numpy() # 疾病相似程度矩阵

dd_association = pd.read_table('./Data/DiDrAMat', sep='\t', header=None).to_numpy() # 药物疾病关联矩阵
dd_association =  np.delete(dd_association, -1, axis=1) # 行疾病列药物
dd_association_t = dd_association.transpose() # 行药物列疾病

alpha = 0.3 # 衰减函数参数
l = 2
r = 2
d = math.log(9999)

drugs_share = np.array(loadmat('./Data/shareWrr.mat')['newWrr']) # 药物和药物之间共同适应疾病的个数
diseases_share = np.array(loadmat('./Data/shareWdd.mat')['newWdd']) # 疾病和疾病之间共同适应药物的个数

c_drug = lib.set_par_fun(dd_association_t, drugs_sim)
c_disease = lib.set_par_fun(dd_association, diseases_sim)

drugs_sim = 1  / (1 + np.exp(c_drug * drugs_sim + d))
diseases_sim = 1  / (1 + np.exp(c_disease * diseases_sim + d))

(drugs_cohesv, diseases_cohesv) = lib.cluster(drugs_sim, diseases_sim, drugs_share, diseases_share)

drugs_cohesv = lib.norm_fun(drugs_cohesv)
diseases_cohesv = lib.norm_fun(diseases_cohesv)

# 随机游走 MBiRW
R0 = dd_association_t / np.concatenate(dd_association_t).sum()
Rt = R0.copy()

for t in range(max(l, r)):
    ftl = 0
    ftr = 0

    if t <= l:
        # m x m * m x n = m x n
        nRtleft = alpha * drugs_cohesv.dot(Rt) + (1 - alpha) * R0
        ftl = 1
    if t <= r:
        # m x n * n * n = m x n
        nRtright = alpha * Rt.dot(diseases_cohesv) + (1 - alpha) * R0
        ftr = 1

    Rt = (ftl * nRtleft + ftr * nRtright) / (ftl + ftr)

drugs_name = pd.read_table('./Data/DrugsName', sep=' ', header=None, squeeze=True).to_numpy()
diseases_name = pd.read_table('./Data/DiseasesName', sep=' ', header=None, squeeze=True).to_numpy()
df = pd.DataFrame(Rt, columns=diseases_name, index=drugs_name)
df.to_csv('./Output/MBiRW.csv')
