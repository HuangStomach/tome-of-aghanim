import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, 
        feature_r1, feature_r2, feature_r3,
        feature_p1, feature_p2, feature_p3):
        super(AutoEncoder, self).__init__()

        # self.encoder_1 = self._gcn(feature_r1[0], feature_r1[1])
        self.encoder_2 = self._gcn(feature_r2[0], feature_r2[1], 0.1)
        self.encoder_3 = self._gcn(feature_r3[0], feature_r3[1], 0.1)
        # self.encoder_4 = self._gcn(feature_r4, feature_r4)

        # self.decoder_1 = self._gcn(feature_p1[0], feature_p1[1])
        self.decoder_2 = self._gcn(feature_p2[0], feature_p2[1], 0.1)
        self.decoder_3 = self._gcn(feature_p3[0], feature_p3[1], 0.1)
        # self.decoder_4 = self._gcn(feature_p4, feature_p4)

        self.encoder = nn.Sequential(
            nn.Linear(feature_r1 + feature_r2[1] + feature_r3[1], 10240),
            nn.Dropout(0.2),
            nn.ReLU(inplace=True),
            # nn.BatchNorm1d(10240),
        )

        self.fc_protein = nn.Linear(10240, feature_p2[0])
        self.fc_disease = nn.Linear(10240, feature_p1)
        
        self.decoder = nn.Sequential(
            nn.Linear(feature_p1 + feature_r1 + feature_p2[1] + feature_p3[1], 10240),
            nn.Dropout(0.2),
            nn.ReLU(inplace=True),
            # nn.BatchNorm1d(10240),
            nn.Linear(10240, feature_r3[0]),
        )

    def _gcn(self, feature_in, feature_out, dropout=0):
        layers = [(GCNConv(feature_in, feature_out), 'x, edge_index -> x1',)]
        if dropout > 0: layers.append(nn.Dropout(dropout))
        layers.append(nn.Sigmoid())
        return Sequential('x, edge_index', layers)

    def forward(self, x1, x2, x3, drug_edge):
        en1 = x1
        en2 = self.encoder_2(x2, drug_edge)
        en3 = self.encoder_3(x3, drug_edge)
        # en4 = self.encoder_4(x4, drug_edge)

        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature)
        rpi_hat = self.fc_protein(encoder0)

        SR_hat = rpi_hat.mm(rpi_hat.t())

        de0 = self.fc_disease(encoder0)
        de1 = x1
        de2 = self.decoder_2(x2, drug_edge)
        de3 = self.decoder_3(x3, drug_edge)
        # de4 = self.encoder_4(x4, drug_edge)

        diease_feature = torch.cat([de0, de1, de2, de3], dim=1)
        decoder0 = self.decoder(diease_feature)

        return rpi_hat, SR_hat, decoder0, decoder0.mm(decoder0.t())
