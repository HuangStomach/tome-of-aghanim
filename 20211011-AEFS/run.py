import logging
import numpy as np
from importlib import import_module
import torch
import torch.nn as nn
import sklearn

import prepare
import metric
from models.autoencoder import AutoEncoder, SONLoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
EPOCH = 100
LR = 0.00001
a1 = 0.0000001
a2 = 0.0000001
formatter = logging.Formatter(
    fmt="%(asctime)s %(levelname)s [%(filename)s] - %(message)s",
    datefmt="%Y-%m-%d %X"
)

def train(model_1, model_2, model_3, dataset):
    # 每次train 要train两组数据 药物和蛋白
    drug_x1 = torch.from_numpy(dataset.drug_x1).float().to(device) # [556, 1024]
    drug_x2 = torch.from_numpy(dataset.drug_x2).float().to(device) # [556, 5603]
    drug_x3 = torch.from_numpy(dataset.drug_x3).float().to(device) # [556, 1512]
    protein_x1 = torch.from_numpy(dataset.protein_x1).float().to(device) # [1512, 128]
    protein_x2 = torch.from_numpy(dataset.protein_x2).float().to(device) # [1512, 5603]
    protein_x3 = torch.from_numpy(dataset.protein_x3).float().to(device) # [1512, 556] 蛋白和药物关联

    SR = dataset.drug_A # 药物相似性
    SP = sklearn.preprocessing.minmax_scale(dataset.protein_A) # 疾病相似性

    [drug_edge, drug_weight] = dataset.edge(dataset.data('rri'), SR)
    drug_edge = torch.from_numpy(drug_edge).long().to(device)
    drug_weight = torch.from_numpy(drug_weight).float().to(device)
    [protein_edge, protein_weight] = dataset.edge(dataset.data('ppi'), SP)
    protein_edge = torch.from_numpy(protein_edge).long().to(device)
    protein_weight = torch.from_numpy(protein_weight).float().to(device)
    
    eye_R = torch.eye(SR.shape[0]).float().to(device)
    eye_P = torch.eye(SP.shape[0]).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)
    SP = torch.from_numpy(SP).float().to(device)

    RPI = torch.from_numpy(dataset.data('rpi')).float().to(device)
    PDI = torch.from_numpy(dataset.data('pdi')).float().to(device)
    
    print("Initialling model...")
    AE = AutoEncoder(
        [1024, 256], [5603, 2048], [1512, 512],
        [128, 128], [5603, 1024], [drug_x1.size()[0], 256],
        drug_num=drug_x1.size()[0], protein_num=protein_x1.size()[0], disease_num=protein_x2.size()[1],
        models=[model_1, model_2, model_3],
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR)
    son_loss = SONLoss(5)
    mse_loss = nn.MSELoss()

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    f = logging.FileHandler('./output/{}_{}_{}.log'.format(model_1, model_2, model_3), encoding='utf-8')
    c = logging.StreamHandler()
    f.setFormatter(formatter)
    c.setFormatter(formatter)
    logger.addHandler(f)
    logger.addHandler(c)

    print("Starting...")
    for epoch in range(EPOCH):
        encoded, SR_hat, decoded, SP_hat = AE(
            drug_x1, drug_x2, drug_x3,
            protein_x1, protein_x2, protein_x3,
            drug_edge, drug_weight, protein_edge, protein_weight
        ) # h3:encoded h6:decoded

        loss1 = mse_loss(encoded, RPI) + son_loss(SR_hat, SR, eye_R, a1)
        loss2 = mse_loss(decoded, PDI) + son_loss(SP_hat, SP, eye_P, a2)

        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        logger.info('[{} {} {}] Epoch: {} train loss: {:.6f}'.format(model_1, model_2, model_3, epoch, loss.item()))

    logger.removeHandler(f)
    logger.removeHandler(c)
    np.savetxt('output/DTINet/{}_{}_{}_RPI.txt'.format(model_1, model_2, model_3), encoded.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/DTINet/{}_{}_{}_PDI.txt'.format(model_1, model_2, model_3), decoded.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.to(device), 'output/DTINet/{}_{}_{}_model.pkl'.format(model_1, model_2, model_3))

if __name__=='__main__':
    while True:
        datasets = ['DTINet', 'AEFS']
        for i, dataset in enumerate(datasets):
            print("[{}] {}".format(i, dataset))
        print("[{}] exit".format(len(datasets)))
        
        str_in = input("Plz select data source: ");
        if str_in.isdigit():
            index = int(str_in)
            if index < 0 or index > len(datasets): continue
            elif index == len(datasets):
                quit()
        else: continue
        dataset = getattr(import_module('dataset'), datasets[index])()
        break

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
            dataset.prepare()
            # models = ['gcn', 'gat', 'gin']
            models = ['gcn', 'gin']

            for model_1 in models:
                for model_2 in models:
                    for model_3 in models:
                        train(model_1, model_2, model_3, dataset)
        elif index == 2:
            metric.run(dataset)
        elif index == 3:
            quit()
        else: continue