import torch
import numpy as np
from sklearn import metrics
import dataset
from models.ae import AutoEncoder

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
def run(type):
    data = dataset.Dataset(type)
    data.init()
    print('model, RPI_auc, RPI_aupr, RPI_aupr_avg, RDI_auc, RDI_aupr, RDI_aupr_avg')
    
    for i in range(10):
        metric(data, i)

def metric(data, tag):
    masks = np.loadtxt('output/{}_masks.txt'.format(tag), dtype=int, delimiter='\n')

    SR = data.drug_A # 药物相似性
    drug_edge, _ = data.edge(SR)
    
    SR = torch.from_numpy(SR).float().to(device)

    AE = AutoEncoder(
        4096, [data.pnum, 1024], [data.dnum, 2048], 
        2048, [data.pnum, 1024], [data.dnum, 2048],
    ).to(device)
    AE_state_dict = torch.load("output/{}_model.pt".format(tag))
    AE.load_state_dict(AE_state_dict)
    AE.eval()

    RPI_hat, _, RDI_hat, _ = AE(
        data.drug_x1, data.drug_x2, data.drug_x3,
        drug_edge
    )

    RPI = data.rpi[masks]
    RDI = data.rdi[masks]
    RPI_hat = RPI_hat.detach().cpu().numpy()[masks]
    RDI_hat = RDI_hat.detach().cpu().numpy()[masks]

    mp = data.metric(RPI, RPI_hat)
    md = data.metric(RDI, RDI_hat)

    line = 'test, {}, {}, {}, {}, {}, {}'.format(
        mp[0], mp[1], mp[2], md[0], md[1], md[2]
    )
    print(line)

    RPI = np.loadtxt('output/{}_RPI.txt'.format(tag))
    RDI = np.loadtxt('output/{}_RDI.txt'.format(tag))
    RPI_hat = np.loadtxt('output/{}_RPI_hat.txt'.format(tag))
    RDI_hat = np.loadtxt('output/{}_RDI_hat.txt'.format(tag))

    mp = data.metric(RPI, RPI_hat)
    md = data.metric(RDI, RDI_hat)

    line = 'train, {}, {}, {}, {}, {}, {}'.format(
        mp[0], mp[1], mp[2], md[0], md[1], md[2]
    )
    print(line)

