# auc aupr
# 我已经跑过了，是85%和4%。我怕自己跑的有问题，所以想让你跑一下验证
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

datasets = ['balance', 'unbalance']
for dataset in datasets:
    label = np.loadtxt('./data/{}/label.out'.format(dataset))
    predict = np.loadtxt('./data/{}/predict.out'.format(dataset))

    fpr, tpr, thresholds = metrics.roc_curve(label, predict, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    print(auc)

    precision, recall, thresholds = metrics.precision_recall_curve(label, predict)
    ap = metrics.average_precision_score(label, predict)
    print(ap)