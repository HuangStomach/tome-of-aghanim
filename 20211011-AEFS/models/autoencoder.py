import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv, GATConv, GINConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, feature_r1, feature_r2, feature_r3, feature_p1, feature_p2, feature_p3, 
        drug_num, protein_num, disease_num, models):
        super(AutoEncoder, self).__init__()
 
        model_1 = getattr(self, '_' + models[0])
        model_2 = getattr(self, '_' + models[1])
        model_3 = getattr(self, '_' + models[2])
        
        self.encoder_1 = model_1(feature_r1[0], feature_r1[1])
        self.encoder_2 = model_2(feature_r2[0], feature_r2[1])
        self.encoder_3 = model_3(feature_r3[0], feature_r3[1])

        self.decoder_1 = model_1(feature_p1[0], feature_p1[1])
        self.decoder_2 = model_2(feature_p2[0], feature_p2[1])
        self.decoder_3 = model_3(feature_p3[0], feature_p3[1])
        
        self.encoder = nn.Sequential(
            nn.Linear(feature_r1[1] + feature_r2[1] + feature_r3[1], 2048),
            # nn.Dropout(0.2),
            nn.Softmax(dim=1),
            nn.BatchNorm1d(2048),
            nn.Linear(2048, protein_num),
            # nn.Dropout(0.2),
            nn.Softmax(dim=1),
        )
        
        self.decoder = nn.Sequential(
            nn.Linear(feature_p1[1] + feature_p2[1] + feature_p3[1], 6144),
            # nn.Dropout(0.2),
            nn.Softmax(dim=1),
            nn.BatchNorm1d(6144),
            nn.Linear(6144, disease_num),
            # nn.Dropout(0.2),
            nn.Softmax(dim=1),
        )
    
    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index', [
            (GCNConv(feature_in, feature_out), 'x, edge_index -> x1'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
            (GCNConv(feature_out, feature_out), 'x1, edge_index -> x2'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
        ])
    
    def _gat(self, feature_in, feature_out, heads=10, dropout=0.2):
        return Sequential('x, edge_index', [
            (GATConv(feature_in, feature_out, heads=heads), 'x, edge_index -> x1'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
            (GATConv(feature_out * heads, feature_out), 'x1, edge_index -> x2'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
        ])
    
    def _gin(self, feature_in, feature_out):
        fc1 = nn.Sequential(
            nn.Linear(feature_in, feature_out),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
            nn.Linear(feature_out, feature_out)
        )

        fc2 = nn.Sequential(
            nn.Linear(feature_out, feature_out),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
            nn.Linear(feature_out, feature_out)
        )

        return Sequential('x, edge_index,', [
            (GINConv(fc1), 'x, edge_index -> x1'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
            (GINConv(fc2), 'x1, edge_index -> x2'),
            # nn.ReLU(inplace=True),
            nn.Softmax(dim=1),
        ])

    def forward(self, r1, r2, r3, p1, p2, p3, drug_edge, drug_weight, protein_edge, protein_weight):
        en1 = self.encoder_1(r1, drug_edge, drug_weight)
        en2 = self.encoder_2(r2, drug_edge, drug_weight)
        en3 = self.encoder_3(r3, drug_edge, drug_weight)
        
        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature) # (batch, protein_num)

        protein_sim = torch.matmul(encoder0.t(), encoder0) # (protein_num, protein_num)
        p1 = protein_sim.matmul(p1)
        p2 = protein_sim.matmul(p2)
        p3 = protein_sim.matmul(p3)

        de1 = self.decoder_1(p1, protein_edge, protein_weight)
        de2 = self.decoder_2(p2, protein_edge, protein_weight)
        de3 = self.decoder_3(p3, protein_edge, protein_weight)

        protein_feature = torch.cat([de1, de2, de3], dim=1)
        decoder0 = self.decoder(protein_feature)

        return encoder0, encoder0.matmul(encoder0.t()), decoder0, decoder0.matmul(decoder0.t())

class SONLoss(nn.Module):
    def __init__(self, k):
        '''
        ## Same order neighbours
        ### Parameters
            * k: top k neighbours
        '''
        super(SONLoss, self).__init__()
        self.k = k
    
    def forward(self, S_hat, S, eye, a):
        '''
        S: 药似性矩阵
        S: 计算出的相似性估计矩阵
        eye: 单位阵
        a: 约束参数
        '''
        _S_hat = S_hat - eye
        _S = S - eye
        S_hat_val, S_hat_idx = _S_hat.topk(self.k)
        S_val, S_idx = _S.topk(self.k)

        diff_idx = S_hat_idx.bitwise_xor(S_idx)
        cal_flag = diff_idx.div(diff_idx).nan_to_num(0) # 索引相同为1 不同为0

        return a * S_hat_val.sub(S_val).pow(2).mul(cal_flag).sum(1).sqrt().sum()