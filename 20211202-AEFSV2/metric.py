import torch
import numpy as np
from sklearn import metrics
import dataset
from models.ae import AutoEncoder

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
def run(dataset, tag='train'):
    file_name = './output/{}_metric.csv'.format(tag)

    dataset.prepare()
    RPI = dataset.rpi
    RDI = dataset.rdi
    print('model, RPI_auc, RPI_aupr, RPI_aupr_avg, RDI_auc, RDI_aupr, RDI_aupr_avg')
    
    with open(file_name, 'w') as f:
        f.write('model, RPI_auc, RPI_aupr, RPI_aupr_avg, RDI_auc, RDI_aupr, RDI_aupr_avg\n')

    with open(file_name, 'a') as f:
        for i in range(5):
            metric(dataset, RPI, RDI, f, i)

def metric(dataset, RPI, RDI, f, tag):
    masks = np.loadtxt('output/{}_masks.txt'.format(tag), dtype=int, delimiter='\n')
    drug_x1 = torch.from_numpy(dataset.drug_x1).float().to(device) # [556, 1024]
    drug_x2 = torch.from_numpy(dataset.drug_x2).float().to(device) # [556, 5603]
    drug_x3 = torch.from_numpy(dataset.drug_x3).float().to(device) # [556, 1512]
    drug_z1 = torch.from_numpy(dataset.drug_z1).float().to(device) # [556, 5603]
    drug_z2 = torch.from_numpy(dataset.drug_z2).float().to(device) # [556, 1512]

    SR = dataset.drug_A # 药物相似性

    drug_edge, drug_weight = dataset.edge(SR, SR)
    drug_edge = torch.from_numpy(drug_edge).long().to(device)
    # drug_weight = torch.from_numpy(drug_weight).float().to(device)
    
    SR = torch.from_numpy(SR).float().to(device)

    AE = AutoEncoder(
        [1024, 1024], [dataset.dnum, 2048], [dataset.pnum, 1024],
        [1024, 1024], [dataset.dnum, 2048], [dataset.pnum, 1024],
        protein_num=dataset.pnum, disease_num=dataset.dnum,
    ).to(device)
    AE_state_dict = torch.load("output/{}_model.pt".format(tag))
    AE.load_state_dict(AE_state_dict)
    AE.eval()

    RPI_hat, _, RDI_hat, _ = AE(
        drug_x1, drug_x2, drug_x3,
        drug_z1, drug_z2, drug_edge
    )

    RPI = dataset.rpi[masks]
    RDI = dataset.rdi[masks]
    RPI_hat = RPI_hat.detach().cpu().numpy()[masks]
    RDI_hat = RDI_hat.detach().cpu().numpy()[masks]

    aupr_rpi_list = []
    for i, row in enumerate(RPI):
        if np.sum(row) == 0 : continue
        aupr_rpi_list.append(metrics.average_precision_score(row, RPI_hat[i]))

    RPI = RPI.flatten()
    RPI_hat = RPI_hat.flatten()

    fpr, tpr, _ = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc_rpi = metrics.auc(fpr, tpr)
    aupr_rpi = metrics.average_precision_score(RPI, RPI_hat)

    aupr_rdi_list = []
    for i, row in enumerate(RDI):
        if np.sum(row) == 0 : continue
        aupr_rdi_list.append(metrics.average_precision_score(row, RDI_hat[i]))

    RDI = RDI.flatten()
    RDI_hat = RDI_hat.flatten()

    fpr, tpr, _ = metrics.roc_curve(RDI, RDI_hat, pos_label=1)
    auc_rdi = metrics.auc(fpr, tpr)
    aupr_rdi = metrics.average_precision_score(RDI, RDI_hat)

    line = 'test, {}, {}, {}, {}, {}, {}'.format(
        auc_rpi, aupr_rpi, np.mean(aupr_rpi_list), 
        auc_rdi, aupr_rdi, np.mean(aupr_rdi_list)
    )
    print(line)
    f.write(line + '\n')

    RPI = np.loadtxt('output/{}_RPI.txt'.format(tag))
    RDI = np.loadtxt('output/{}_RDI.txt'.format(tag))
    RPI_hat = np.loadtxt('output/{}_RPI_hat.txt'.format(tag))
    RDI_hat = np.loadtxt('output/{}_RDI_hat.txt'.format(tag))

    aupr_rpi_list = []
    for i, row in enumerate(RPI):
        if np.sum(row) == 0 : continue
        aupr_rpi_list.append(metrics.average_precision_score(row, RPI_hat[i]))

    RPI = RPI.flatten()
    RPI_hat = RPI_hat.flatten()

    fpr, tpr, _ = metrics.roc_curve(RPI, RPI_hat, pos_label=1)
    auc_rpi = metrics.auc(fpr, tpr)
    aupr_rpi = metrics.average_precision_score(RPI, RPI_hat)

    aupr_rdi_list = []
    for i, row in enumerate(RDI):
        if np.sum(row) == 0 : continue
        aupr_rdi_list.append(metrics.average_precision_score(row, RDI_hat[i]))

    RDI = RDI.flatten()
    RDI_hat = RDI_hat.flatten()

    fpr, tpr, _ = metrics.roc_curve(RDI, RDI_hat, pos_label=1)
    auc_rdi = metrics.auc(fpr, tpr)
    aupr_rdi = metrics.average_precision_score(RDI, RDI_hat)

    line = 'train, {}, {}, {}, {}, {}, {}'.format(
        auc_rpi, aupr_rpi, np.mean(aupr_rpi_list), 
        auc_rdi, aupr_rdi, np.mean(aupr_rdi_list)
    )
    print(line)
    f.write(line + '\n')

if __name__=='__main__':
    dataset = dataset.Dataset()
    run(dataset)
