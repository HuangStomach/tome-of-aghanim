import numpy as np
from sklearn import metrics
import dataset

file_name = 'metric.csv'
def run(dataset):
    dataset.prepare()
    RPI = dataset.rpi
    
    print('model, RPI_auc, RPI_aupr(RPI_aupr_avg)')
    with open(file_name, 'w') as f:
        f.write('model, RPI_auc, RPI_aupr(RPI_aupr_avg)\n')

    with open(file_name, 'a') as f:
        metric(RPI, f)

def metric(RPI, f, tag='train'): 
    RPI_hat = np.loadtxt('output/DTINet/{}_RPI.txt'.format(tag))

    aupr_rpi_list = []
    for i, row in enumerate(RPI):
        if np.sum(row) == 0 : continue
        aupr_rpi_list.append(metrics.average_precision_score(row, RPI_hat[i]))

    RPI = RPI.flatten()
    RPI_hat = RPI_hat.flatten()

    fpr, tpr, _ = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc_rpi = metrics.auc(fpr, tpr)
    aupr_rpi = metrics.average_precision_score(RPI, RPI_hat)

    print('{}, {}({})'.format(auc_rpi, aupr_rpi, np.mean(aupr_rpi_list)))
    f.write('{}, {}({})\n'.format(auc_rpi, aupr_rpi, np.mean(aupr_rpi_list)))

if __name__=='__main__':
    dataset = dataset.Dataset()
    dataset.prepare()
    RPI = dataset.rpi
    print('model, RPI_auc, RPI_aupr(RPI_aupr_avg)')
    
    with open(file_name, 'w') as f:
        f.write('model, RPI_auc, RPI_aupr(RPI_aupr_avg)\n')

    with open(file_name, 'a') as f:
        for i in range(5):
            metric(RPI, f, i)