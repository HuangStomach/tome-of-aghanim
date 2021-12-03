import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv, GATConv, GINConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, feature_r1, feature_r2, feature_r3, feature_p1, feature_p2, feature_p3, protein_num, disease_num):
        super(AutoEncoder, self).__init__()
        
        self.encoder_1 = self._gcn(feature_r1[0], feature_r1[1])
        self.encoder_2 = self._gcn(feature_r2[0], feature_r2[1])
        self.encoder_3 = self._gcn(feature_r3[0], feature_r3[1])

        self.decoder_1 = self._gcn(feature_p1[0], feature_p1[1])
        self.decoder_2 = self._gcn(feature_p2[0], feature_p2[1])
        self.decoder_3 = self._gcn(feature_p3[0], feature_p3[1])
        
        self.encoder = nn.Sequential(
            nn.Linear(feature_r1[1] + feature_r2[1] + feature_r3[1], 2048),
            # nn.Dropout(0.2),
            nn.Relu(),
            nn.Linear(2048, protein_num),
            # nn.Dropout(0.2),
            # nn.Tanh(),
            # nn.BatchNorm1d(protein_num)
        )
        
        self.decoder = nn.Sequential(
            nn.Linear(feature_p1[1] + feature_p2[1] + feature_p3[1], 6144),
            # nn.Dropout(0.2),
            nn.ReLU(dim=1),
            nn.Linear(6144, disease_num),
            # nn.Dropout(0.2),
            nn.Sigmoid(dim=1),
        )
    
    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index, edge_weight', [
            (GCNConv(feature_in, feature_out), 'x, edge_index, edge_weight -> x1'),
            nn.ReLU(inplace=True),
            # nn.Softmax(dim=1),
            (GCNConv(feature_out, feature_out), 'x1, edge_index, edge_weight -> x2'),
            nn.ReLU(inplace=True),
            # nn.Softmax(dim=1),
        ])

    def forward(self, r1, r2, r3, p1, p2, SP, drug_edge, drug_weight, protein_edge, protein_weight):
        en1 = self.encoder_1(r1, drug_edge, edge_weight=drug_weight)
        en2 = self.encoder_2(r2, drug_edge, edge_weight=drug_weight)
        en3 = self.encoder_3(r3, drug_edge, edge_weight=drug_weight)

        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature) # (batch, protein_num)

        protein_sim = torch.matmul(encoder0.t(), encoder0) # (protein_num, protein_num)
        p3 = SP.matmul(encoder0.t())

        de1 = self.decoder_1(p1, protein_edge, edge_weight=protein_weight)
        de2 = self.decoder_2(p2, protein_edge, edge_weight=protein_weight)
        de3 = self.decoder_3(p3, protein_edge, edge_weight=protein_weight)

        protein_feature = torch.cat([de1, de2, de3], dim=1)
        decoder0 = self.decoder(protein_feature)

        return encoder0, protein_sim, decoder0, decoder0.matmul(decoder0.t())
