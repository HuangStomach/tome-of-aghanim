from importlib import import_module
import torch
from torch.utils.data import TensorDataset, DataLoader
from lib import *
import protein

device = 'cuda' if torch.cuda.is_available() else 'cpu'
EPOCH = 200
BATCH_SIZE = 64
LR = 0.0001
# drug_num = 1307
# protein_num = 1996
# indication_num = 3926
drug_feature = 1024  # ECFPs指纹
a1 = 0.00000001
a2 = 0.0001

def train(_model, dataset):
    # 每次train 要train两组数据 药物和蛋白
    Model = getattr(import_module('models.{}'.format(_model.lower())), _model)()

    print("初始化模型")
    drug_set = TensorDataset(dataset.id, dataset.drug_x1, dataset.drug_x2, dataset.drug_x3)
    protein_set = TensorDataset(dataset.protein_x1, dataset.protein_x2, dataset.protein_x3)
    drug_loader = DataLoader(dataset=drug_set, batch_size=BATCH_SIZE, shuffle=True, num_workers=6)
    protein_loader = DataLoader(dataset=protein_set, batch_size=BATCH_SIZE, shuffle=True, num_workers=6)

    drug_model = Model(dataset.drug_x1.size() , )
    AE = AutoEncoder(drug_feature, indication_num)
    optimizer = torch.optim.Adam(AE.parameters(), lr=LR)
    aefs_loss = SONLoss(SR, 10)
    mse_loss = nn.MSELoss()
    print("cuda加速")
    AE = AE.to(device)
    SP = SP.to(device)
    print("开始训练")
    for epoch in range(EPOCH):
        for i, data in enumerate(train_loader):
            batch_id, batch_x, batch_h, batch_y = data
            batch_SR = torch.empty(batch_x.shape[0], batch_x.shape[0])
            for m in range(batch_x.shape[0]):
                for n in range(batch_x.shape[0]):
                    batch_SR[m, n] = SR[batch_id[m], batch_id[n]]
            batch_SR = batch_SR.to(device)
            batch_x = batch_x.to(device)
            batch_h = batch_h.to(device)
            batch_y = batch_y.to(device)

            encoded, decoded = AE(batch_x) # h3:encoded h6:decoded
            loss1 = mse_loss(encoded, batch_h) + aefs_loss(encoded, a1)
            loss2 = mse_loss(decoded, batch_y) + aefs_loss(decoded, a2)
            loss = loss1 + loss2
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        print('Epoch:', epoch, 'train loss: %.20f' % loss.to(device).data)
    torch.save(AE.to(device), 'result/model.pkl')

while True:
    datasets = ['DTINet', 'AEFS']
    for i, dataset in enumerate(datasets):
        print("[{}] {}".format(i, dataset))
    
    str_in = input("请选择数据源: ");
    if str_in.isdigit():
        index = int(str_in)
        if index < 0 or index >= len(datasets): continue
    else: continue
    dataset = getattr(import_module('dataset'), datasets[index])()
    break

while True:
    print("[0] 处理蛋白序列")
    print("[1] 随机切分数据")
    print("[2] 训练数据")
    print("[3] 退出")
    str_in = input("请选择需要进行的操作: ");
    if not str_in.isdigit(): continue

    index = int(str_in)

    if index == 0:
        protein.transform()
        print('\033[32m完成\033[0m')
    elif index == 1:
        dataset.split_data()
        print('\033[32m完成\033[0m')
    elif index == 2:
        dataset.prepare()
        models = ['GAT', 'GCN', 'GIN']

        for model_1 in models:
            for model_2 in models:
                for model_3 in models:
                    train(model_1, dataset)
                    train(model_2, dataset)
                    train(model_3, dataset)
    elif index == 3:
        quit()
    else: continue