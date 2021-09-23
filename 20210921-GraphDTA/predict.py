import torch
import numpy as np
from utils import *
from models.ginconv import GINConvNet

TEST_BATCH_SIZE = 512
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def predicting(model, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            output = model(data)
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat(
                (total_labels, data.y.view(-1, 1).cpu()), 0)
    return total_labels.numpy().flatten(), total_preds.numpy().flatten()


datasets = ['balance', 'unbalance']
epochs = [1000, 100]
for i, dataset in enumerate(datasets):
    modeling = torch.load('./data/{}/model_GIN_epoch_{}.pkl'.format(dataset, epochs[i]))
    model_st = str(modeling)

    test = TestbedDataset(root='data', dataset=dataset + '_test')
    test_loader = DataLoader(test, batch_size=TEST_BATCH_SIZE, shuffle=False)
    model = modeling.to(device)
    G, P = predicting(model, test_loader)
    np.savetxt('./data/{}/label.out'.format(dataset), G)
    np.savetxt('./data/{}/predict.out'.format(dataset), P)