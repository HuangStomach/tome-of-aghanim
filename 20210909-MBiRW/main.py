import math
import pandas as pd
import numpy as np
from scipy.io import loadmat 
import lib

drugs_sim = pd.read_table('./Data/DrugSimMat', sep=' ', header=None).to_numpy() # 药物相似程度矩阵
drugs_name = pd.read_table('./Data/DrugsName', sep=' ', header=None, squeeze=True).to_numpy()

diseases_sim = pd.read_table('./Data/DiseaseSimMat', sep=' ', header=None).to_numpy() # 疾病相似程度矩阵
diseases_name = pd.read_table('./Data/DiseasesName', sep=' ', header=None, squeeze=True).to_numpy()

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

(a, b) = lib.cluster(drugs_sim, diseases_sim, drugs_share, diseases_share, drugs_name, diseases_name)