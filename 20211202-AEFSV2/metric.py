import torch
import numpy as np
from sklearn import metrics
import dataset
from models.ae import AutoEncoder

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
def run(dataset):
    dataset = dataset.Dataset()
    dataset.prepare()
    print('model, RPI_auc, RPI_aupr, RPI_aupr_avg, RDI_auc, RDI_aupr, RDI_aupr_avg')
    
    for i in range(10):
        metric(dataset, i)

def metric(dataset, tag):
    masks = np.loadtxt('output/{}_masks.txt'.format(tag), dtype=int, delimiter='\n')

    SR = dataset.drug_A # 药物相似性
    drug_edge, _ = dataset.edge(SR)
    
    SR = torch.from_numpy(SR).float().to(device)

    AE = AutoEncoder(
        4096, [dataset.dnum, 2048], [dataset.pnum, 1024], [4192, 128],
        2048, [dataset.dnum, 2048], [dataset.pnum, 1024], [4192, 128],
    ).to(device)
    AE_state_dict = torch.load("output/{}_model.pt".format(tag))
    AE.load_state_dict(AE_state_dict)
    AE.eval()

    RPI_hat, _, RDI_hat, _ = AE(
        dataset.drug_x1, dataset.drug_x2, dataset.drug_x3, dataset.drug_x4,
        drug_edge
    )

    RPI = dataset.rpi[masks]
    RDI = dataset.rdi[masks]
    RPI_hat = RPI_hat.detach().cpu().numpy()[masks]
    RDI_hat = RDI_hat.detach().cpu().numpy()[masks]

    mp = dataset.metric(RPI, RPI_hat)
    md = dataset.metric(RDI, RDI_hat)

    line = 'test, {}, {}, {}, {}, {}, {}'.format(
        mp[0], mp[1], mp[2], md[0], md[1], md[2]
    )
    print(line)

    RPI = np.loadtxt('output/{}_RPI.txt'.format(tag))
    RDI = np.loadtxt('output/{}_RDI.txt'.format(tag))
    RPI_hat = np.loadtxt('output/{}_RPI_hat.txt'.format(tag))
    RDI_hat = np.loadtxt('output/{}_RDI_hat.txt'.format(tag))

    mp = dataset.metric(RPI, RPI_hat)
    md = dataset.metric(RDI, RDI_hat)

    line = 'train, {}, {}, {}, {}, {}, {}'.format(
        mp[0], mp[1], mp[2], md[0], md[1], md[2]
    )
    print(line)

if __name__=='__main__':
    dataset = dataset.Dataset()
    run(dataset)
