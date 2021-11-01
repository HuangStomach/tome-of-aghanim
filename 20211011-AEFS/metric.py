import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

def run(dataset):
    models = ['gcn', 'gat', 'gin']
    RPI = dataset.data('rpi')
    PDI = dataset.data('pdi')

    for model_1 in models:
        for model_2 in models:
            for model_3 in models:
                metric(model_1, model_2, model_3, RPI, PDI)

def metric(model_1, model_2, model_3, RPI, PDI): 
    RPI_hat = np.loadtxt('output/DTINet/{}_{}_{}_RPI.txt'.format(model_1, model_2, model_3))
    PDI_hat = np.loadtxt('output/DTINet/{}_{}_{}_PDI.txt'.format(model_1, model_2, model_3))

    fpr, tpr, thresholds = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    print('{}_{}_{}_{}_auc: '.format(model_1, model_2, model_3, 'RPI'), auc)

    aupr = metrics.average_precision_score(RPI, RPI_hat)
    print('{}_{}_{}_{}_aupr: '.format(model_1, model_2, model_3, 'RPI'), aupr)

    fpr, tpr, thresholds = metrics.roc_curve(PDI, PDI_hat, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    print('{}_{}_{}_{}_auc: '.format(model_1, model_2, model_3, 'PDI'), auc)

    aupr = metrics.average_precision_score(PDI, PDI_hat)
    print('{}_{}_{}_{}_aupr: '.format(model_1, model_2, model_3, 'PDI'), aupr)
