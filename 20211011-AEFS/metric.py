import numpy as np
from sklearn import metrics
from lib import *

file_name = 'metric.csv'
def run(dataset):
    # models = ['gcn', 'gat', 'gin']
    models = ['gcn']
    dataset.prepare()
    RPI = dataset.rpi
    PDI = dataset.pdi
    
    print('model, RPI_auc, RPI_aupr, PDI_auc, PDI_aupr')
    with open(file_name, 'w') as f:
        f.write('model, RPI_auc, RPI_aupr, PDI_auc, PDI_aupr\n')

    with open(file_name, 'a') as f:
        for model_1 in models:
            for model_2 in models:
                for model_3 in models:
                    metric(model_1, model_2, model_3, RPI, PDI, f)

def metric(model_1, model_2, model_3, RPI, PDI, f): 
    RPI_hat = np.loadtxt('output/DTINet/{}_{}_{}_RPI.txt'.format(model_1, model_2, model_3))
    PDI_hat = np.loadtxt('output/DTINet/{}_{}_{}_PDI.txt'.format(model_1, model_2, model_3))

    # along_drugs = []
    # for i in range(RPI.shape[0]):
    #     if np.sum(RPI[i]) == 0 and np.sum(RPI_hat[i]): along_drugs.append(i)

    # along_proteins = []
    # for j in range(RPI.shape[1]):
    #     if np.sum(RPI[:,j]) == 0 and np.sum(RPI_hat[:, j]): along_proteins.append(j)

    # RPI = np.delete(RPI, along_drugs, axis=0)
    # RPI = np.delete(RPI, along_proteins, axis=1)
    # RPI_hat = np.delete(RPI_hat, along_drugs, axis=0)
    # RPI_hat = np.delete(RPI_hat, along_proteins, axis=1)

    count = 0
    aupr_rpi_list = []
    for i, row in enumerate(RPI):
        if np.sum(row) == 0 : continue

        aupr_rpi_list.append(metrics.average_precision_score(row, RPI_hat[i]))
    
    print(np.average(aupr_rpi_list))

    RPI = RPI.flatten()
    RPI_hat = RPI_hat.flatten()

    PDI = PDI.flatten()
    PDI_hat = PDI_hat.flatten()

    fpr, tpr, thresholds = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc_rpi = metrics.auc(fpr, tpr)
    aupr_rpi = metrics.average_precision_score(RPI, RPI_hat)

    fpr, tpr, thresholds = metrics.roc_curve(PDI, PDI_hat, pos_label=1)
    auc_pdi = metrics.auc(fpr, tpr)
    aupr_pdi = metrics.average_precision_score(PDI, PDI_hat)

    print('{}_{}_{}, {}, {}, {}, {}'.format(model_1, model_2, model_3, auc_rpi, aupr_rpi, auc_pdi, aupr_pdi))
    f.write('{}_{}_{}, {}, {}, {}, {}\n'.format(model_1, model_2, model_3, auc_rpi, aupr_rpi, auc_pdi, aupr_pdi))

# def metric(model_1, model_2, model_3, RPI, PDI, f): 
#     RPI_hat = np.loadtxt('output/DTINet/{}_{}_{}_RPI.txt'.format(model_1, model_2, model_3))
#     PDI_hat = np.loadtxt('output/DTINet/{}_{}_{}_PDI.txt'.format(model_1, model_2, model_3))
#     print(metrics.roc_auc_score(RPI, RPI_hat))

#     auc_list = []
#     aupr_list = []
#     tpr_list = []
#     fpr_list = []
#     recall_list = []
#     precision_list = []
#     c = 0
#     for i in range(len(RPI)):
#         if np.sum(RPI[i]) == 0:
#             c += 1
#             continue
#         else:
#             tpr1, fpr1, precision1, recall1 = tpr_fpr_precision_recall(RPI[i], RPI_hat[i])
#             fpr_list.append(fpr1)
#             tpr_list.append(tpr1)
#             precision_list.append(precision1)
#             recall_list.append(recall1)
#             auc_list.append(metrics.auc(fpr1, tpr1))
#             aupr_list.append(metrics.auc(recall1, precision1)+recall1[0]*precision1[0])

#     coverage = []
#     for i in tpr_list:
#         try:
#             coverage.append(i.index(1.0)+1)
#         except:
#             print('1')

#     tpr = equal_len_list(tpr_list)
#     fpr = equal_len_list(fpr_list)
#     precision = equal_len_list(precision_list)
#     recall = equal_len_list(recall_list)
#     tpr_mean = np.mean(tpr, axis=0)
#     fpr_mean = np.mean(fpr, axis=0)
#     recall_mean = np.mean(recall, axis=0)
#     precision_mean = np.mean(precision, axis=0)

#     print('{}_{}_{}, {}, {}'.format(model_1, model_2, model_3, metrics.roc_auc_score(RPI, RPI_hat), metrics.auc(recall_mean, precision_mean)+recall_mean[0]*precision_mean[0]))
