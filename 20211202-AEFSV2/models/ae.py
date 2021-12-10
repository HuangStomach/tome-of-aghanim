import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, feature_r1, feature_r2, feature_r3, feature_p1, feature_p2, feature_p3, protein_num, disease_num):
        super(AutoEncoder, self).__init__()
        self.protein_num = protein_num
        
        self.encoder_1 = self._gcn(feature_r1[0], feature_r1[1])
        self.encoder_2 = self._gcn(feature_r2[0], feature_r2[1])
        self.encoder_3 = self._gcn(feature_r3[0], feature_r3[1])

        self.decoder_1 = self._gcn(feature_p1[0], feature_p1[1])
        self.decoder_2 = self._gcn(feature_p2[0], feature_p2[1])
        self.decoder_3 = self._gcn(feature_p3[0], feature_p3[1])
        
        self.encoder = nn.Sequential(
            nn.Linear(feature_r1[1] + feature_r2[1] + feature_r3[1], 8192),
            nn.Dropout(0.1),
            nn.LeakyReLU(inplace=True),
            # nn.Tanh(),
            nn.Linear(8192, protein_num + 1024),
        )
        
        self.decoder = nn.Sequential(
            nn.Linear(feature_p1[1] + feature_p2[1] + feature_p3[1], 8192),
            nn.Dropout(0.1),
            nn.LeakyReLU(inplace=True),
            # nn.Tanh(),
            nn.Linear(8192, disease_num),
            # nn.Sigmoid(dim=1),
        )
    
    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index, edge_weight', [
            (GCNConv(feature_in, feature_out), 'x, edge_index, edge_weight -> x1',),
            nn.Dropout(0.1),
            # nn.LeakyReLU(inplace=True),
            nn.Sigmoid(),
            # nn.Softmax(dim=1),
            (GCNConv(feature_out, feature_out), 'x1, edge_index, edge_weight -> x2'),
            nn.Dropout(0.1),
            nn.Sigmoid(),
            # nn.LeakyReLU(inplace=True),
            # nn.Softmax(dim=1),
        ])

    def forward(self, x1, x2, x3, z1, z2, SR, drug_edge, drug_weight):
        en1 = self.encoder_1(x1, drug_edge, edge_weight=drug_weight)
        en2 = self.encoder_2(x2, drug_edge, edge_weight=drug_weight)
        en3 = self.encoder_3(x3, drug_edge, edge_weight=drug_weight)

        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature) # (batch, protein_num)
        rpi_hat = encoder0[:, :self.protein_num]
        hidden_state = encoder0[:, self.protein_num:]

        drug_sim_1 = rpi_hat.matmul(rpi_hat.t()) # (protein_num, protein_num)

        de1 = self.decoder_1(SR.matmul(hidden_state), drug_edge, edge_weight=drug_weight)
        de2 = self.decoder_2(z1, drug_edge, edge_weight=drug_weight)
        de3 = self.decoder_3(z2, drug_edge, edge_weight=drug_weight)

        diease_feature = torch.cat([de1, de2, de3], dim=1)
        decoder0 = self.decoder(diease_feature)

        return rpi_hat, drug_sim_1, decoder0, decoder0.matmul(decoder0.t())

