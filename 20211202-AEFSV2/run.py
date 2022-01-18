import time
import numpy as np
import torch
import torch.nn as nn
from sklearn import metrics
import logging
from logging import handlers

import metric
import dataset
from models.ae import AutoEncoder
from models.loss import SONLoss, WeightMSELoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
EPOCH = 2000
LR = 0.001  # 0.00009
WD = 0.000004 # 0.000004

sim_threshold = 0.5
loss_p_weight = 0.998
loss_d_weight = 0.95
# loss_p_weight = 0.995
# loss_d_weight = 0.95
a1 = 0.00000001
a2 = 0.00000001

def train(trainData, testData, mask, logger, tag='train'):
    # 每次train 要train两组数据 药物和蛋白
    drug_x1 = torch.from_numpy(trainData.drug_x1).float().to(device) # [556, 1024]
    drug_x2 = torch.from_numpy(trainData.drug_x2).float().to(device) # [556, 5603]
    drug_x3 = torch.from_numpy(trainData.drug_x3).float().to(device) # [556, 1512]
    drug_x4 = torch.from_numpy(trainData.drug_x4).float().to(device) # [556, 1512]
    drug_z1 = torch.from_numpy(trainData.drug_z1).float().to(device) # [556, 5603]
    drug_z2 = torch.from_numpy(trainData.drug_z2).float().to(device) # [556, 1512]

    drug_x1_test = torch.from_numpy(testData.drug_x1).float().to(device) # [556, 1024]
    drug_x2_test = torch.from_numpy(testData.drug_x2).float().to(device) # [556, 5603]
    drug_x3_test = torch.from_numpy(testData.drug_x3).float().to(device) # [556, 1512]
    drug_x4_test = torch.from_numpy(testData.drug_x4).float().to(device) # [556, 1512]
    drug_z1_test = torch.from_numpy(testData.drug_z1).float().to(device) # [556, 5603]
    drug_z2_test = torch.from_numpy(testData.drug_z2).float().to(device) # [556, 1512]

    SR = trainData.drug_A # 药物相似性
    SR_test = testData.drug_A # 药物相似性

    drug_edge, drug_weight = trainData.edge(SR, sim_threshold)
    drug_edge = torch.from_numpy(drug_edge).long().to(device)
    # drug_weight = torch.from_numpy(drug_weight).float().to(device)
    drug_edge_test, drug_weight_test = testData.edge(SR_test, sim_threshold)
    drug_edge_test = torch.from_numpy(drug_edge_test).long().to(device)
    # drug_weight = torch.from_numpy(drug_weight).float().to(device)

    eye_R = torch.eye(trainData.rnum).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)

    RPI = torch.from_numpy(trainData.rpi).float().to(device)
    RDI = torch.from_numpy(trainData.rdi).float().to(device)
    RPI_test = testData.rpi[mask]
    RDI_test = testData.rdi[mask]
    RPI_test_f = RPI_test.flatten()

    print("Initialling model...")
    AE = AutoEncoder(
        [1024, 1024], [trainData.dnum, 2048], [trainData.pnum, 1024], [4192, 128],
        [2048, 2048], [trainData.dnum, 2048], [trainData.pnum, 1024],
        protein_num=trainData.pnum, disease_num=trainData.dnum,
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR, weight_decay=WD)
    son_loss = SONLoss(10)
    mse_loss_p = WeightMSELoss(loss_p_weight)
    mse_loss_d = WeightMSELoss(loss_d_weight)
    # mse_loss_p = nn.MSELoss()
    # mse_loss_d = nn.MSELoss()

    print("Starting {}...".format(tag))
    for epoch in range(EPOCH):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(
            drug_x1, drug_x2, drug_x3, drug_x4,
            drug_z1, drug_z2, drug_edge
        ) # h3:encoded h6:decoded

        loss1 = mse_loss_p(RPI_hat, RPI)# + son_loss(SR_hat_1, SR, eye_R, a1)
        loss2 = mse_loss_d(RDI_hat, RDI)# + son_loss(SR_hat_2, SR, eye_R, a2)

        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        RPI_hat_test, _, RDI_hat_test, _ = AE(
            drug_x1_test, drug_x2_test, drug_x3_test, drug_x4_test,
            drug_z1_test, drug_z2_test, drug_edge_test
        ) # h3:encoded h6:decoded
        
        RPI_hat_test = RPI_hat_test.detach().cpu().numpy()[mask]

        aupr_rpi_list = []
        for i, row in enumerate(RPI_test):
            if np.sum(row) == 0 : continue
            aupr_rpi_list.append(metrics.average_precision_score(row, RPI_hat_test[i]))

        RPI_hat_test = RPI_hat_test.flatten()
        fpr, tpr, _ = metrics.roc_curve(RPI_test_f, RPI_hat_test, pos_label=1)
        auc_rpi = metrics.auc(fpr, tpr)
        aupr_rpi = metrics.average_precision_score(RPI_test_f, RPI_hat_test)

        info = 'Epoch: {} train loss: {:.6f}, test auc: {:.6f}, aupr: {:.6f}, aupr_mean: {:.6f}'.format(
            epoch, loss.item(), auc_rpi, aupr_rpi, np.mean(aupr_rpi_list)
        )
        logger.info(info)
        # print()

    np.savetxt('output/{}_RPI_hat.txt'.format(tag), RPI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RPI.txt'.format(tag), RPI.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI_hat.txt'.format(tag), RDI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI.txt'.format(tag), RDI.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.state_dict(), 'output/{}_model.pt'.format(tag))

if __name__=='__main__':
    trainData = dataset.Dataset()
    filename = './output/{}.log'.format(
        time.strftime("%Y%m%d%H%M%S", time.localtime())
    )
    logger = logging.getLogger(filename)
    logger.setLevel(logging.DEBUG)

    sh = logging.StreamHandler()
    logger.addHandler(sh)
    fh = handlers.RotatingFileHandler(filename=filename)
    logger.addHandler(fh)

    while True:
        print("[0] train")
        print("[1] metric")
        print("[2] exit")
        str_in = input("Plz select the opt: ");
        if not str_in.isdigit(): continue

        index = int(str_in)
        if index == 0:
            drug_count = trainData.drugs().shape[0]
            shuffled_drugs = np.arange(drug_count)
            drugs = np.arange(drug_count)
            np.random.shuffle(shuffled_drugs)
            splits = np.array_split(shuffled_drugs, 10)
            
            for i in range(10):
                trainData.prepare(mask_drugs=splits[i])

                testData = dataset.Dataset()
                testData.prepare()

                train(trainData, testData, splits[i], logger, i)
                np.savetxt('output/{}_masks.txt'.format(i), splits[i], fmt='%d')
        elif index == 1:
            metric.run(trainData)
        elif index == 2:
            quit()
        else: continue
