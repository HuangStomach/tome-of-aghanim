import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, 
        feature_r1, feature_r2, feature_r3, feature_r4,
        feature_p1, feature_p2, feature_p3, 
        protein_num, disease_num):

        super(AutoEncoder, self).__init__()
        self.protein_num = protein_num

        # self.encoder_1 = self._gcn(feature_r1[0], feature_r1[1])
        self.encoder_2 = self._gcn(feature_r2[0], feature_r2[1])
        self.encoder_3 = self._gcn(feature_r3[0], feature_r3[1])
        self.encoder_4 = self._gcn(feature_r4[0], feature_r4[1])

        # self.decoder_1 = self._gcn(feature_p1[0], feature_p1[1])
        self.decoder_2 = self._gcn(feature_p2[0], feature_p2[1])
        self.decoder_3 = self._gcn(feature_p3[0], feature_p3[1])

        self.encoder = nn.Sequential(
            nn.Linear(feature_r1[1] + feature_r2[1] + feature_r3[1] + feature_r4[1], 8196),
            nn.ReLU(inplace=True),
            nn.Dropout(0.4),
        )

        self.fc_protein = nn.Linear(8196, protein_num)
        self.fc_disease = nn.Linear(8196, feature_p1[1])

        self.decoder = nn.Sequential(
            nn.Linear(feature_p1[1] + feature_p2[1] + feature_p3[1], 10240),
            nn.ReLU(inplace=True),
            nn.Dropout(0.4),
            nn.Linear(10240, disease_num),
        )

    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index', [
            (GCNConv(feature_in, feature_out), 'x, edge_index -> x1',),
            nn.Sigmoid(),
        ])

    def forward(self, x1, x2, x3, x4, z1, z2, drug_edge):
        en1 = x1
        en2 = self.encoder_2(x2, drug_edge)
        en3 = self.encoder_3(x3, drug_edge)
        en4 = self.encoder_4(x4, drug_edge)

        drug_feature = torch.cat([en1, en2, en3, en4], dim=1)
        encoder0 = self.encoder(drug_feature)
        rpi_hat = self.fc_protein(encoder0)

        SR_hat = rpi_hat.matmul(rpi_hat.t())

        de1 = self.fc_disease(encoder0)
        de2 = self.decoder_2(z1, drug_edge)
        de3 = self.decoder_3(z2, drug_edge)
        # de4 = x1

        diease_feature = torch.cat([de1, de2, de3], dim=1)
        decoder0 = self.decoder(diease_feature)

        return rpi_hat, SR_hat, decoder0, decoder0.matmul(decoder0.t())
