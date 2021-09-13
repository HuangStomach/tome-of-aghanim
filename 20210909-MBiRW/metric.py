import math
import pandas as pd
import numpy as np
from scipy.io import loadmat 
np.set_printoptions(suppress = True)

dd_association = pd.read_table('./Data/DiDrAMat', sep='\t', header=None).to_numpy() # 药物疾病关联矩阵
dd_association =  np.delete(dd_association, -1, axis=1).transpose() # 行疾病列药物
dd_association_pre = pd.read_csv('./Output/MBiRW.csv', header=None).to_numpy()

(r, c) = dd_association_pre.shape
flatten = dd_association_pre.flatten()
flatten.sort()
percentages = np.array([1, 10000, 90000, 185400])
marks = np.zeros((4, 2))

for k in range(percentages.shape[0]):
    threshold = flatten[percentages[k]]
    TP = 0
    FN = 0
    FP = 0
    TN = 0
    for i in range(r):
        for j in range(c):
            value = 1 if dd_association_pre[i][j] > threshold else 0
            if value == 1 and dd_association[i][j] == 1: TP += 1
            if value == 0 and dd_association[i][j] == 0: TN += 1
            if value == 0 and dd_association[i][j] == 1: FN += 1
            if value == 1 and dd_association[i][j] == 0: FP += 1

    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    marks[k] = [TPR, FPR]

auc = 0.0
for i in range(marks.shape[0] - 1):
    auc += (marks[i + 1][0] - marks[i][0]) * (marks[i][1] + marks[i + 1][1])
print(auc)