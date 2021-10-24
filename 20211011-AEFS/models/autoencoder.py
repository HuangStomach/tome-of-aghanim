import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import Sequential, GCNConv, GATConv, GINConv, global_max_pool as gmp
import numpy as np
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, feature_d1, feature_d2, feature_d3, feature_p1, feature_p2, feature_p3, 
        drug_num, protein_num, disease_num, models):
        super(AutoEncoder, self).__init__()
        
        model_1 = getattr(self, '_' + models[0])
        model_2 = getattr(self, '_' + models[1])
        model_3 = getattr(self, '_' + models[2])
        
        self.encoder_1 = model_1(feature_d1[0], feature_d1[1])
        self.encoder_2 = model_2(feature_d2[0], feature_d2[1])
        self.encoder_3 = model_3(feature_d3[0], feature_d3[1])

        self.decoder_1 = model_1(feature_p1[0], feature_p1[1])
        self.decoder_2 = model_2(feature_p2[0], feature_p2[1])
        self.decoder_3 = model_3(feature_p3[0], feature_p3[1])
        
        self.encoder = nn.Sequential(
            nn.Linear(feature_d1[1] + feature_d2[1] + feature_d3[1], 2048),
            nn.Dropout(0.2),
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
            nn.Linear(feature_p1[1] + feature_p2[1] + feature_p3[1], 6144),
            nn.Dropout(0.2),
            nn.GELU(),
            nn.BatchNorm1d(6144),
            nn.Linear(6144, 6144),
            nn.Dropout(0.2),
            nn.GELU(),
            nn.BatchNorm1d(6144),
            nn.Linear(6144, disease_num),
            nn.Dropout(0.2),
        )
    
    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index', [
            (GCNConv(feature_in, feature_in), 'x, edge_index -> x1'),
            nn.ReLU(),
            (nn.Dropout(0.2), 'x1 -> x1'),
            (GCNConv(feature_in, feature_out), 'x1, edge_index -> x2'),
            nn.ReLU(),
            (nn.Dropout(0.2), 'x2 -> x2'),
        ])
    
    def _gat(self, feature_in, feature_out, heads=10, dropout=.2):
        return Sequential('x, edge_index', [
            (GATConv(feature_in, feature_in, heads=heads, dropout=dropout), 'x, edge_index -> x1'),
            nn.ReLU(),
            (GATConv(feature_in, feature_out, dropout=dropout), 'x1, edge_index -> x2'),
            nn.ReLU(),
        ])
    
    def _gin(self, feature_in, feature_out):
        fc = nn.Sequential(
            nn.Linear(feature_in, feature_out),
            nn.ReLU(),
            nn.Linear(feature_out, feature_out)
        )

        return Sequential('x, edge_index', [
            (GINConv(fc), 'x, edge_index -> x1'),
            nn.ReLU(),
            (GINConv(fc), 'x1, edge_index -> x2'),
            nn.ReLU(),
        ])

    def forward(self, d1, d2, d3, p1, p2, p3, drug_edge, protein_edge):
        en1 = self.encoder_1(d1, drug_edge)
        en2 = self.encoder_2(d2, drug_edge)
        en3 = self.encoder_3(d3, drug_edge)
        
        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature) # (batch, protein_num)
        encoder1 = F.softmax(encoder0, dim=1)

        protein_sim = torch.matmul(encoder1.t(), encoder1) # (protein_num, protein_num)
        p1 = protein_sim.matmul(p1, protein_edge)
        p2 = protein_sim.matmul(p2, protein_edge)
        p3 = protein_sim.matmul(p3, protein_edge)

        de1 = self.decoder_1(p1)
        de2 = self.decoder_2(p2)
        de3 = self.decoder_3(p3)

        protein_feature = torch.cat([de1, de2, de3], dim=1)
        decoder0 = self.decoder(protein_feature)
        decoder1 = F.softmax(decoder0, dim=1)
        return encoder1, decoder1

class SONLoss(nn.Module):
    def __init__(self, k):
        super(SONLoss, self).__init__()
        self.k = k
    
    def forward(self, S, S_hat, eye, a):
        '''
        S: 药似性矩阵
        S: 计算出的相似性估计矩阵
        eye: 单位阵
        a: 约束参数
        '''
        S_hat_val, S_hat_idx = S_hat.topk(self.k) - eye
        S_val, S_idx = S.topk(self.k) - eye

        diff_idx = S_hat_idx.bitwise_xor(S_idx)
        punish = diff_idx.abs().sum(1) # 对位近邻的索引之差视为惩罚项
        cal_flag = diff_idx.div(diff_idx).nan_to_num(0) # 索引相同为1 不同为0
        
        return a * S_hat_val.sub(S_val).pow(2).mul(cal_flag).sum(1).sqrt().mul(punish).sum()