import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, feature_r1, feature_r2, feature_r3, feature_p1, feature_p2, feature_p3, protein_num, disease_num):
        super(AutoEncoder, self).__init__()
        self.protein_num = protein_num

        # self.encoder_1 = self._gcn(feature_r1[0], feature_r1[1])
        self.encoder_2 = self._gcn(feature_r2[0], feature_r2[1])
        self.encoder_3 = self._gcn(feature_r3[0], feature_r3[1])

        # self.decoder_1 = self._gcn(feature_p1[0], feature_p1[1])
        self.decoder_2 = self._gcn(feature_p2[0], feature_p2[1])
        self.decoder_3 = self._gcn(feature_p3[0], feature_p3[1])

        self.encoder = nn.Sequential(
            nn.Linear(feature_r1[1] + feature_r2[1] + feature_r3[1], 8192),
            nn.LeakyReLU(inplace=True),
            nn.Dropout(0.6),
            # nn.Linear(8192, protein_num + 1024),
        )
        self.fc_p = nn.Linear(8192, protein_num)
        self.fc_d = nn.Linear(8192, 2048)

        self.decoder = nn.Sequential(
            nn.Linear(2048 + feature_p2[1] + feature_r1[1], 10240),
            nn.LeakyReLU(inplace=True),
            nn.Dropout(0.6),
            nn.Linear(10240, disease_num),
        )

    def _gcn(self, feature_in, feature_out):
        return Sequential('x, edge_index', [
            (GCNConv(feature_in, feature_out), 'x, edge_index -> x1',),
            nn.Sigmoid(),
            nn.Dropout(0.5),
            # (GCNConv(feature_out, feature_out), 'x1, edge_index, edge_weight -> x2'),
            # nn.Sigmoid(),
        ])

    def forward(self, x1, x2, x3, z1, z2, SR, drug_edge, drug_weight):
        # en1 = self.encoder_1(x1, drug_edge, edge_weight=drug_weight)
        en1 = x1
        en2 = self.encoder_2(x2, drug_edge)
        en3 = self.encoder_3(x3, drug_edge)

        drug_feature = torch.cat([en1, en2, en3], dim=1)
        encoder0 = self.encoder(drug_feature)
        rpi_hat = self.fc_p(encoder0)
        hidden_state = self.fc_d(encoder0)

        SR_hat = rpi_hat.matmul(rpi_hat.t())

        # de1 = self.decoder_1(SR.matmul(hidden_state), drug_edge, edge_weight=drug_weight)
        de1 = hidden_state
        de2 = self.decoder_2(z1, drug_edge)
        #de3 = self.decoder_3(z2, drug_edge)
        de4 = x1

        diease_feature = torch.cat([de1, de2, de4], dim=1)
        decoder0 = self.decoder(diease_feature)

        return rpi_hat, SR_hat, decoder0, decoder0.matmul(decoder0.t())
