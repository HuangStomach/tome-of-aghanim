import torch
import torch.nn as nn
from torch_geometric.nn import Sequential, GCNConv
from lib import *

class AutoEncoder(nn.Module):
    def __init__(self, params):
        super(AutoEncoder, self).__init__()
        self.fr = fr = params['fr_dim'] # [ecfps, rpi, rdi, rgo]
        self.fd = fd = params['fd_dim'] # [ecfps, rpi, rdi, rgo]

        self.encoder_2 = self._gcn(fr[1][0], fr[1][1], params['graph_dropout'])
        self.encoder_3 = self._gcn(fr[2][0], fr[2][1], params['graph_dropout'])
        if len(fr) > 3: self.encoder_4 = self._gcn(fr[3][0], fr[3][1], params['graph_dropout'])

        self.decoder_2 = self._gcn(fd[1][0], fd[1][1], params['graph_dropout'])
        self.decoder_3 = self._gcn(fd[2][0], fd[2][1], params['graph_dropout'])
        if len(fd) > 3: self.encoder_4 = self._gcn(fd[3][0], fd[3][1], params['graph_dropout'])

        f_in = fr[0] + fr[1][1] + fr[2][1] 
        if len(fr) > 3: f_in += fr[3][1]
        self.encoder = nn.Sequential(
            nn.Linear(f_in, 10240),
            nn.Dropout(params['dropout']),
            nn.ReLU(inplace=True),
        )

        self.fc_protein = nn.Linear(10240, fd[1][0])
        self.fc_disease = nn.Linear(10240, fd[0])

        f_in = fd[0] + fr[0] + fd[1][1] + fd[2][1]
        if len(fd) > 3: f_in += fd[3][1]
        self.decoder = nn.Sequential(
            nn.Linear(f_in, 10240),
            nn.Dropout(params['dropout']),
            nn.ReLU(inplace=True),
            nn.Linear(10240, fr[2][0]),
        )

    def _gcn(self, f_in, f_out, dropout=0.0):
        layers = [(GCNConv(f_in, f_out), 'x, edge_index -> x1',)]
        if dropout > 0: layers.append(nn.Dropout(dropout))
        layers.append(nn.Sigmoid())
        return Sequential('x, edge_index', layers)

    def forward(self, data):
        en = [data.drug_x1]
        en.append(self.encoder_2(data.drug_x2, data.drug_edge))
        en.append(self.encoder_3(data.drug_x3, data.drug_edge))
        if len(self.fr) > 3: en.append(self.encoder_4(data.drug_x4, data.drug_edge))

        drug_feature = torch.cat(en, dim=1)
        encoder0 = self.encoder(drug_feature)
        rpi_hat = self.fc_protein(encoder0)

        SR_hat = rpi_hat.mm(rpi_hat.t())

        de = [self.fc_disease(encoder0)]
        de.append(data.drug_x1)
        de.append(self.decoder_2(data.drug_x2, data.drug_edge))
        de.append(self.decoder_3(data.drug_x3, data.drug_edge))
        if len(self.fd) > 3: de.append(self.decoder_4(data.x4, data.drug_edge))

        diease_feature = torch.cat(de, dim=1)
        decoder0 = self.decoder(diease_feature)

        return rpi_hat, SR_hat, decoder0, decoder0.mm(decoder0.t())
