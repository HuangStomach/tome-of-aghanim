import time
import numpy as np
import torch
import torch.nn as nn
import logging
from logging import handlers

import metric
import dataset
from models.ae import AutoEncoder
from models.loss import SONLoss, WeightMSELoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
EPOCH = 1000
LR = 0.00003  # 0.00009
WD = 0.000004 # 0.000004

sim_threshold = 0.5
loss_p_weight = 0.998
loss_d_weight = 0.95
# loss_p_weight = 0.995
# loss_d_weight = 0.95
a1 = 0.00000001
a2 = 0.00000001

def train(trainData, testData, mask, logger, tag='train'):
    SR = trainData.drug_A # 药物相似性
    SR_test = testData.drug_A # 药物相似性

    drug_edge, drug_weight = trainData.edge(SR, sim_threshold)
    drug_edge_test, drug_weight_test = testData.edge(SR_test, sim_threshold)

    eye_R = torch.eye(trainData.rnum).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)

    RPI = torch.from_numpy(trainData.rpi).float().to(device)
    RDI = torch.from_numpy(trainData.rdi).float().to(device)
    RPI_test = testData.rpi[mask]
    RDI_test = testData.rdi[mask]

    print("Initialling model...")
    AE = AutoEncoder(
        4096, [trainData.dnum, 2048], [trainData.pnum, 1024], [4192, 128],
        2048, [trainData.dnum, 2048], [trainData.pnum, 1024], [4192, 128],
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR, weight_decay=WD)
    son_loss = SONLoss(8)
    mse_loss_p = WeightMSELoss(loss_p_weight)
    mse_loss_d = WeightMSELoss(loss_d_weight)
    # mse_loss_p = nn.MSELoss()
    # mse_loss_d = nn.MSELoss()

    print("Starting {}...".format(tag))
    for epoch in range(EPOCH):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(
            trainData.drug_x1, trainData.drug_x2, trainData.drug_x3, trainData.drug_x4,
            drug_edge
        )

        loss1 = mse_loss_p(RPI_hat, RPI) + son_loss(SR_hat_1, SR, eye_R, a1)
        loss2 = mse_loss_d(RDI_hat, RDI) + son_loss(SR_hat_2, SR, eye_R, a2)

        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        RPI_hat_test, _, RDI_hat_test, _ = AE(
            testData.drug_x1, testData.drug_x2, testData.drug_x3, testData.drug_x4,
            drug_edge_test
        )
        
        RPI_hat_test = RPI_hat_test.detach().cpu().numpy()[mask]
        RDI_hat_test = RDI_hat_test.detach().cpu().numpy()[mask]
        mp = testData.metric(RPI_test, RPI_hat_test)
        md = testData.metric(RDI_test, RDI_hat_test)

        info = 'Epoch: {} loss: {:.6f}, pauc: {:.6f}, paupr: {:.6f}, paupr_m: {:.6f}, dauc: {:.6f}, daupr: {:.6f}, daupr_m: {:.6f}'.format(
            epoch, loss.item(), mp[0], mp[1], mp[2], md[0], md[1], md[2]
        )
        logger.info(info)

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
