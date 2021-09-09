import math
import pandas as pd
import numpy as np
from scipy.io import loadmat 
import lib

drugs_sim = pd.read_table('./Data/DrugSimMat', sep=' ', header=None) # 药物相似程度矩阵
drugs_name = pd.read_table('./Data/DrugsName', sep=' ', header=None).to_numpy()[:, 0]

diseases_sim = pd.read_table('./Data/DiseaseSimMat', sep=' ', header=None) # 疾病相似程度矩阵
diseases_name = pd.read_table('./Data/DiseasesName', sep=' ', header=None).to_numpy()[:, 0]

dd_association = pd.read_table('./Data/DiDrAMat', sep=' ', header=None) # 药物疾病关联矩阵
dd_association_t = dd_association.transpose()

alpha = 0.3 # 衰减函数参数
l = 2
r = 2
d = math.log(9999)

drugs_share = np.array(loadmat('./Data/shareWrr.mat')['newWrr'])
diseases_share = np.array(loadmat('./Data/shareWdd.mat')['newWdd'])
lib.set_par_fun(dd_association_t, drugs_sim)
