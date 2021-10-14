import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from torch_geometric.nn import GCNConv, global_max_pool as gmp
import numpy as np
from lib import *

device = 'cuda' if torch.cuda.is_available() else 'cpu'
EPOCH = 200
BATCH_SIZE = 64
LR = 0.0001
drug_num = 1307
protein_num = 1996
indication_num = 3926
drug_feature = 1024  # ECFPs指纹
a1 = 0.00000001
a2 = 0.0001

class GCN(nn.Module):
    def __init__(self, size_x, size_y):
        super(GCN, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(size_x, 2048),
            nn.Dropout(0.2),         # 参数是扔掉的比例
            nn.GELU(),
            nn.BatchNorm1d(2048),
            nn.Linear(2048, 2048),
            nn.Dropout(0.2),
            nn.GELU(),
            nn.BatchNorm1d(2048),
            nn.Linear(2048, protein_num),
            nn.Dropout(0.2),
        )

        self.decoder = nn.Sequential(
            nn.Linear(protein_num, 4096),
            nn.Dropout(0.2),
            nn.GELU(),
            nn.BatchNorm1d(4096),
            nn.Linear(4096, 4096),
            nn.Dropout(0.2),
            nn.GELU(),
            nn.BatchNorm1d(4096),
            nn.Linear(4096, size_y),
            nn.Dropout(0.2),
        )

    # def forward(self, x, sp):
    def forward(self, x):
        e0 = self.encoder(x)
        e1 = F.softmax(e0, dim=1)
        d0 = self.decoder(e1)
        d1 = F.softmax(d0, dim=1)
        return e1, d1

