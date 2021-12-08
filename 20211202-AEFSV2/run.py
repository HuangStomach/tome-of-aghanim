import numpy as np
import torch
import torch.nn as nn
import dataset

import prepare
import metric
from models.ae import AutoEncoder
from models.loss import SONLoss, FocalLoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
EPOCH = 1000
LR = 0.00001
a1 = 0.00000001
a2 = 0.00000001

def train(dataset, tag='train'):
    # 每次train 要train两组数据 药物和蛋白
    drug_x1 = torch.from_numpy(dataset.drug_x1).float().to(device) # [556, 1024]
    drug_x2 = torch.from_numpy(dataset.drug_x2).float().to(device) # [556, 5603]
    drug_x3 = torch.from_numpy(dataset.drug_x3).float().to(device) # [556, 1512]
    drug_z1 = torch.from_numpy(dataset.drug_z1).float().to(device) # [556, 5603]
    drug_z2 = torch.from_numpy(dataset.drug_z2).float().to(device) # [556, 1512]

    SR = dataset.drug_A # 药物相似性

    drug_edge, drug_weight = dataset.edge(SR, SR)
    drug_edge = torch.from_numpy(drug_edge).long().to(device)
    drug_weight = torch.from_numpy(drug_weight).float().to(device)
    
    eye_R = torch.eye(dataset.rnum).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)

    RPI = torch.from_numpy(dataset.rpi).float().to(device)
    RDI = torch.from_numpy(dataset.rdi).float().to(device)
    
    print("Initialling model...")
    AE = AutoEncoder(
        [1024, 512], [dataset.dnum, 1024], [dataset.pnum, 512],
        [1024, 1024], [dataset.dnum, 1024], [dataset.pnum, 512],
        protein_num=dataset.pnum, disease_num=dataset.dnum,
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR)
    son_loss = SONLoss(10)
    mse_loss = nn.MSELoss()

    print("Starting...")
    for epoch in range(EPOCH):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(
            drug_x1, drug_x2, drug_x3,
            drug_z1, drug_z2, SR,
            drug_edge, drug_weight
        ) # h3:encoded h6:decoded

        loss1 = mse_loss(RPI_hat, RPI) + son_loss(SR_hat_1, SR, eye_R, a1)
        loss2 = mse_loss(RDI_hat, RDI) + son_loss(SR_hat_2, SR, eye_R, a2)

        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        print('Epoch: {} train loss: {:.6f}'.format(epoch, loss.item()))

    np.savetxt('output/{}_RPI.txt'.format(tag), RPI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI.txt'.format(tag), RDI_hat.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.to(device), 'output/{}_model.pkl'.format(tag))

if __name__=='__main__':
    dataset = dataset.Dataset()

    while True:
        print("[0] protein embedding")
        print("[1] train")
        print("[2] metric")
        print("[3] exit")
        str_in = input("Plz select the opt: ");
        if not str_in.isdigit(): continue

        index = int(str_in)

        if index == 0:
            prepare.proteins()
            print('\033[32mfinish.\033[0m')
        elif index == 1:
            drugs = np.arange(dataset.data('rri').shape[0])
            np.random.shuffle(drugs)
            splits = np.array_split(drugs, 5)
            for i in range(5):
                dataset.prepare(mask_drugs=splits[i])
                train(dataset, i)
        elif index == 2:
            metric.run(dataset)
        elif index == 3:
            quit()
        else: continue