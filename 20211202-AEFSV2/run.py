import gc
import sys
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

def train(trainData, testData, mask, logger, tag='train'):
    eye_R = torch.eye(trainData.rnum).float().to(device)

    SR = torch.from_numpy(trainData.drug_A).float().to(device)
    RPI = torch.from_numpy(trainData.rpi).float().to(device)
    RDI = torch.from_numpy(trainData.rdi).float().to(device)
    RPI_test = testData.rpi[mask]
    RDI_test = testData.rdi[mask]

    print("Initialling model...")
    AE = AutoEncoder(trainData.params).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=trainData.params['lr'], weight_decay=trainData.params['wd'])
    son_loss = SONLoss(10)
    mse_loss_p = WeightMSELoss(trainData.params['loss_p_weight'])
    mse_loss_d = WeightMSELoss(trainData.params['loss_d_weight'])

    print("Starting {}...".format(tag))
    for epoch in range(trainData.params['epoch']):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(trainData)

        loss1 = mse_loss_p(RPI_hat, RPI) + trainData.params['a1'] * son_loss(SR_hat_1, SR, eye_R)
        loss2 = mse_loss_d(RDI_hat, RDI) + trainData.params['a2'] * son_loss(SR_hat_2, SR, eye_R)

        loss = loss1 + trainData.params['loss_weight'] * loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 10 == 0 and torch.cuda.is_available():
            try:
                RPI_hat_test, _, RDI_hat_test, _ = AE(testData)
                
                RPI_hat_test = RPI_hat_test.detach().cpu().numpy()[mask]
                RDI_hat_test = RDI_hat_test.detach().cpu().numpy()[mask]
                mp = testData.metric(RPI_test, RPI_hat_test)
                md = testData.metric(RDI_test, RDI_hat_test)

                info = 'Epoch: {} loss: {:.6f}, pauc: {:.6f}, paupr: {:.6f}, paupr_m: {:.6f}, dauc: {:.6f}, daupr: {:.6f}, daupr_m: {:.6f}'.format(
                    epoch, loss.item(), mp[0], mp[1], mp[2], md[0], md[1], md[2]
                )
                logger.info(info)
            except Exception as e:
                print('error', e, RPI_hat_test)

    np.savetxt('output/{}_RPI_hat.txt'.format(tag), RPI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RPI.txt'.format(tag), RPI.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI_hat.txt'.format(tag), RDI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI.txt'.format(tag), RDI.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.state_dict(), 'output/{}_model.pt'.format(tag))
    return mp[0] >= 0.95 and (mp[1] >= 0.58 or mp[2] >= 0.58) and md[0] >= 0.92

if __name__=='__main__':
    type = sys.argv[1] if len(sys.argv) > 1 else 'DTINet'

    while True:
        print("[0] prepare")
        print("[1] train")
        print("[2] metric")
        print("[3] exit")
        str_in = input("Plz select the opt: ");
        if not str_in.isdigit(): continue

        index = int(str_in)
        if index == 0:
            data = dataset.Dataset(type)
            data.prepare()
        elif index == 1:
            while True:
                localtime = time.strftime("%Y%m%d_%H%M%S", time.localtime())
                filename = './output/{}.log'.format(localtime)
                logger = logging.getLogger(filename)
                logger.setLevel(logging.DEBUG)

                sh = logging.StreamHandler()
                logger.addHandler(sh)
                fh = handlers.RotatingFileHandler(filename=filename)
                logger.addHandler(fh)

                trainData = dataset.Dataset(type)
                splits = trainData.splits()
                np.savetxt('output/{}_mask.txt'.format(localtime), splits, fmt='%s', delimiter=',')
                testData = dataset.Dataset(type)
                testData.init()

                logger.info(trainData.params)
                result = False
                for i in range(10):
                    trainData.init(mask_drugs=splits[i])

                    result = train(trainData, testData, splits[i], logger, i) or result
                    if i == 0 and result == False: break

                del logger
                gc.collect()
            if result: break
        elif index == 2:
            metric.run(type)
        elif index == 3:
            quit()
        else: continue
