import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

file_name = 'metric.csv'
def run(dataset):
    models = ['gcn', 'gat', 'gin']
    RPI = dataset.data('rpi').flatten()
    PDI = dataset.data('pdi').flatten()
    print('model, RPI_auc, RPI_aupr, PDI_auc, PDI_aupr')
    with open(file_name, 'w') as f:
        f.write('model, RPI_auc, RPI_aupr, PDI_auc, PDI_aupr\n')

    with open(file_name, 'a') as f:
        for model_1 in models:
            for model_2 in models:
                for model_3 in models:
                    metric(model_1, model_2, model_3, RPI, PDI, f)

def metric(model_1, model_2, model_3, RPI, PDI, f): 
    RPI_hat = np.loadtxt('output/DTINet/{}_{}_{}_RPI.txt'.format(model_1, model_2, model_3)).flatten()
    PDI_hat = np.loadtxt('output/DTINet/{}_{}_{}_PDI.txt'.format(model_1, model_2, model_3)).flatten()
 
    fpr, tpr, thresholds = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc_rpi = metrics.auc(fpr, tpr)
    aupr_rpi = metrics.average_precision_score(RPI, RPI_hat)

    fpr, tpr, thresholds = metrics.roc_curve(PDI, PDI_hat, pos_label=1)
    auc_pdi = metrics.auc(fpr, tpr)
    aupr_pdi = metrics.average_precision_score(PDI, PDI_hat)

    print('{}_{}_{}, {}, {}, {}, {}'.format(model_1, model_2, model_3, auc_rpi, aupr_rpi, auc_pdi, aupr_pdi))
    f.write('{}_{}_{}, {}, {}, {}, {}\n'.format(model_1, model_2, model_3, auc_rpi, aupr_rpi, auc_pdi, aupr_pdi))
