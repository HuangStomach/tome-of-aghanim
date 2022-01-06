import numpy as np
import torch
import torch.nn as nn
import dataset

import metric
from models.ae import AutoEncoder
from models.loss import SONLoss, WeightMSELoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
EPOCH = 1000
LR = 0.0001 # 0.00009
WD = 0.000004 # 0.000004

loss_p_weight = 0.999
loss_d_weight = 0.99
a1 = 0.00000001
a2 = 0.00000001

def train(trainData, testData, tag='train'):
    # 每次train 要train两组数据 药物和蛋白
    drug_x1 = torch.from_numpy(trainData.drug_x1).float().to(device) # [556, 1024]
    drug_x2 = torch.from_numpy(trainData.drug_x2).float().to(device) # [556, 5603]
    drug_x3 = torch.from_numpy(trainData.drug_x3).float().to(device) # [556, 1512]
    drug_z1 = torch.from_numpy(trainData.drug_z1).float().to(device) # [556, 5603]
    drug_z2 = torch.from_numpy(trainData.drug_z2).float().to(device) # [556, 1512]

    drug_x1_test = torch.from_numpy(testData.drug_x1).float().to(device) # [556, 1024]
    drug_x2_test = torch.from_numpy(testData.drug_x2).float().to(device) # [556, 5603]
    drug_x3_test = torch.from_numpy(testData.drug_x3).float().to(device) # [556, 1512]
    drug_z1_test = torch.from_numpy(testData.drug_z1).float().to(device) # [556, 5603]
    drug_z2_test = torch.from_numpy(testData.drug_z2).float().to(device) # [556, 1512]

    SR = trainData.drug_A # 药物相似性
    SR_test = testData.drug_A # 药物相似性

    drug_edge, drug_weight = trainData.edge(SR, SR)
    drug_edge = torch.from_numpy(drug_edge).long().to(device)
    # drug_weight = torch.from_numpy(drug_weight).float().to(device)
    drug_edge_test, drug_weight_test = testData.edge(SR_test, SR_test)
    drug_edge_test = torch.from_numpy(drug_edge_test).long().to(device)
    # drug_weight = torch.from_numpy(drug_weight).float().to(device)

    eye_R = torch.eye(trainData.rnum).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)

    RPI = torch.from_numpy(trainData.rpi).float().to(device)
    RDI = torch.from_numpy(trainData.rdi).float().to(device)
    RPI_test = torch.from_numpy(testData.rpi).float().to(device)
    RDI_test = torch.from_numpy(testData.rdi).float().to(device)

    print("Initialling model...")
    AE = AutoEncoder(
        [1024, 1024], [trainData.dnum, 2048], [trainData.pnum, 1024],
        [1024, 1024], [trainData.dnum, 2048], [trainData.pnum, 1024],
        protein_num=trainData.pnum, disease_num=trainData.dnum,
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR, weight_decay=WD)
    son_loss = SONLoss(10)
    mse_loss_p = WeightMSELoss(loss_p_weight)
    mse_loss_d = WeightMSELoss(loss_d_weight)

    mse_loss_p_test = nn.MSELoss()
    mse_loss_d_test = nn.MSELoss()

    print("Starting {}...".format(tag))
    for epoch in range(EPOCH):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(
            drug_x1, drug_x2, drug_x3,
            drug_z1, drug_z2, drug_edge
        ) # h3:encoded h6:decoded

        loss1 = mse_loss_p(RPI_hat, RPI) + son_loss(SR_hat_1, SR, eye_R, a1)
        loss2 = mse_loss_d(RDI_hat, RDI) + son_loss(SR_hat_2, SR, eye_R, a2)

        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        RPI_hat_test, _, RDI_hat_test, _ = AE(
            drug_x1_test, drug_x2_test, drug_x3_test,
            drug_z1_test, drug_z2_test, drug_edge_test
        ) # h3:encoded h6:decoded
        loss_test = mse_loss_p_test(RPI_hat_test, RPI_test) + mse_loss_d_test(RDI_hat_test, RDI_test)

        print('Epoch: {} train loss: {:.6f}, test loss: {:.6f}'.format(epoch, loss.item(), loss_test.item()))

    np.savetxt('output/{}_RPI_hat.txt'.format(tag), RPI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RPI.txt'.format(tag), RPI.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI_hat.txt'.format(tag), RDI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI.txt'.format(tag), RDI.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.state_dict(), 'output/{}_model.pt'.format(tag))

if __name__=='__main__':
    trainData = dataset.Dataset()

    while True:
        print("[0] train")
        print("[1] metric")
        print("[2] exit")
        str_in = input("Plz select the opt: ");
        if not str_in.isdigit(): continue

        index = int(str_in)
        if index == 0:
            drugs = np.arange(trainData.data('rri').shape[0])
            np.random.shuffle(drugs)
            splits = np.array_split(drugs, 5)
            for i in range(5):
                trainData.prepare(mask_drugs=splits[i])

                testData = dataset.Dataset()
                testData.prepare(mask_drugs=np.delete(drugs, splits[i], axis=0))

                train(trainData, testData, i)
                np.savetxt('output/{}_masks.txt'.format(i), splits[i], fmt='%d')
        elif index == 1:
            metric.run(trainData)
        elif index == 2:
            quit()
        else: continue
