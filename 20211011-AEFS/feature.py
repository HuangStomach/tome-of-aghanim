import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

device = 'cuda' if torch.cuda.is_available() else 'cpu'
BATCH_SIZE = 64
LR = 0.0005
NUM_EPOCHS = 1000

def output(model_d, model_p, dataset):
    drug_A = dataset.drug_sim('train')
    drug_x1 = drug_A * dataset.fps('train')
    drug_x2 = drug_A * dataset.dpi('train')
    drug_x3 = drug_A * dataset.rda('train')

    protein_A = dataset.protein_sim('train')
    protein_x1 = protein_A * dataset.fps('train').T
    protein_x2 = protein_A * dataset.dpi('train').T
    protein_x3 = protein_A * dataset.rda('train').T
    # model
    model_drug = model_d().to(device)
    model_protein = model_p().to(device)
    # optimizer
    optimizer_drug = torch.optim.Adam(model_drug.parameters(), lr=LR)
    optimizer_protein = torch.optim.Adam(model_protein.parameters(), lr=LR)
    # data
    train_set = TensorDataset(drug_x1)
    train_set = TensorDataset(drug_x2)
    train_set = TensorDataset(drug_x3)

    train_set = TensorDataset(protein_x1)
    train_set = TensorDataset(protein_x2)
    train_set = TensorDataset(protein_x3)

    train_loader = DataLoader(dataset=train_set, batch_size=BATCH_SIZE, shuffle=True, num_workers=6)

    for epoch in range(NUM_EPOCHS):
        output = train(model_drug, device, train_loader, optimizer_drug, epoch + 1)
        output = train(model_protein, device, train_loader, optimizer_protein, epoch + 1)

        print('Epoch:{} saving output'.format(epoch))
        

def train(model, device, data_loader, optimizer, epoch):
    for data in data_loader:
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)

        loss_fn = nn.MSELoss()
        loss = loss_fn(output, data.y.view(-1, 1).float().to(device))
        loss.backward()
        optimizer.step()

        print('Epoch:', epoch, 'train loss: %.20f' % loss.to(device).data)

    return output