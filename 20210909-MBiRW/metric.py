import math
import pandas as pd
import numpy as np
from scipy.io import loadmat 
np.set_printoptions(suppress = True)

association = np.array(loadmat('./Data/Datasets_indep/DiDrMat.mat')['R_Wdr'])
drugs_name = pd.read_table('./Data/Datasets_indep/R_Wrname', sep=' ', header=None, squeeze=True).to_numpy()
diseases_name = pd.read_table('./Data/Datasets_indep/R_Wdname', sep=' ', header=None, squeeze=True).to_numpy()
association = pd.DataFrame(association.transpose(), columns=diseases_name, index=drugs_name)
prediction = pd.read_csv('./Output/MBiRW.csv', index_col=0)

(r, c) = association.shape
flatten = prediction.to_numpy().flatten()[::-1]
flatten.sort()
percentages = np.array([0.01, 0.1, 0.2, 0.5, 0.8, 0.95, 1])
marks = np.zeros((percentages.shape[0], 2))

for k in range(percentages.shape[0]):
    threshold = flatten[min(flatten.shape[0] - 1, int(flatten.shape[0] * percentages[k]))]
    TP = FN = FP = TN = 0
    for i in association.index:
        for j in association.columns:
            if i not in prediction.index or j not in prediction.columns: continue
            value = 1 if prediction.loc[i, j] > threshold else 0
            if value == 1 and association.loc[i, j] == 1: TP += 1
            if value == 0 and association.loc[i, j] == 0: TN += 1
            if value == 0 and association.loc[i, j] == 1: FN += 1
            if value == 1 and association.loc[i, j] == 0: FP += 1

    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    marks[k] = [TPR, FPR]

auc = 0.0
marks = marks[::-1]

for i in range(marks.shape[0] - 1):
    auc += (marks[i + 1][0] - marks[i][0]) * (marks[i][1] + marks[i + 1][1])
