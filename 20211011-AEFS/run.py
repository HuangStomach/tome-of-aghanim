from importlib import import_module
from itertools import cycle
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from lib import *
import prepare
from models.autoencoder import AutoEncoder, SONLoss

device = 'cuda' if torch.cuda.is_available() else 'cpu'
EPOCH = 200
BATCH_SIZE = 4
LR = 0.0001
# drug_num = 1307
# protein_num = 1996
# indication_num = 3926
drug_feature = 1024  # ECFPs指纹
a1 = 0.00000001
a2 = 0.0001

def train(model_1, model_2, model_3, dataset):
    # 每次train 要train两组数据 药物和蛋白
    drug_x1 = torch.from_numpy(dataset.drug_x1[:10, :]).float() # [556, 1024]
    drug_x2 = torch.from_numpy(dataset.drug_x2[:10, :]).float() # [556, 5603]
    drug_x3 = torch.from_numpy(dataset.drug_x3[:10, :]).float() # [556, 1512]
    protein_x1 = torch.from_numpy(dataset.protein_x1[:20, :]).float() # [1512, 128]
    protein_x2 = torch.from_numpy(dataset.protein_x2[:20, :]).float() # [1512, 5603]
    protein_x3 = torch.from_numpy(dataset.protein_x3[:20, :]).float() # [1512, 556]

    SR = dataset.drug_A[:10, :10] # 药物相似性
    SP = max_min_normalize(dataset.protein_sim()[:20, :20]) # 疾病相似性
    drug_edge = torch.from_numpy(dataset.edge_index(SR)).long()
    protein_edge = torch.from_numpy(dataset.edge_index(SP)).long()
    
    eye_R = torch.eye(SR.shape[0]).float()
    eye_P = torch.eye(SP.shape[0]).float()
    SR = torch.from_numpy(SR).float()
    SP = torch.from_numpy(SP).float()

    print("初始化模型")
    AE = AutoEncoder(
        [1024, 128], [5603, 1024], [1512, 256],
        [128, 128], [1512, 1024], [556, 256],
        drug_num=drug_x1.size()[0], protein_num=protein_x1.size()[0], disease_num=protein_x2.size()[1],
        models=[model_1, model_2, model_3],
    )
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR)
    son_loss = SONLoss(5)
    mse_loss = nn.MSELoss()

    print("开始训练")
    for epoch in range(EPOCH):
        encoded, decoded = AE(
            drug_x1, drug_x2, drug_x3,
            protein_x1, protein_x2, protein_x3,
            drug_edge, protein_edge
        ) # h3:encoded h6:decoded

        loss1 = mse_loss(encoded, SR) + son_loss(SR, encoded, eye_R, a1)
        loss2 = mse_loss(decoded, SP) + son_loss(SP, decoded, eye_P, a2)
        loss = loss1 + loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        print('Epoch:', epoch, 'train loss: %.20f' % loss.to(device).data)

    torch.save(AE.to(device), 'output/DTINet/model.pkl')

if __name__=='__main__':
    while True:
        datasets = ['DTINet', 'AEFS']
        for i, dataset in enumerate(datasets):
            print("[{}] {}".format(i, dataset))
        print("[{}] 退出".format(len(datasets)))
        
        str_in = input("请选择数据源: ");
        if str_in.isdigit():
            index = int(str_in)
            if index < 0 or index > len(datasets): continue
            elif index == len(datasets):
                quit()
        else: continue
        dataset = getattr(import_module('dataset'), datasets[index])()
        break

    while True:
        print("[0] 处理蛋白序列")
        print("[1] 生成疾病相似性")
        print("[2] 随机切分数据")
        print("[3] 训练数据")
        print("[4] 退出")
        str_in = input("请选择需要进行的操作: ");
        if not str_in.isdigit(): continue

        index = int(str_in)

        if index == 0:
            prepare.proteins()
            print('\033[32m完成\033[0m')
        if index == 1:
            prepare.diseases(dataset)
            print('\033[32m完成\033[0m')
        elif index == 2:
            dataset.split_data()
            print('\033[32m完成\033[0m')
        elif index == 3:
            dataset.prepare()
            models = ['gat', 'gcn', 'gin']

            for model_1 in models:
                for model_2 in models:
                    for model_3 in models:
                        train(model_1, model_2, model_3, dataset)
                        quit()
        elif index == 4:
            quit()
        else: continue