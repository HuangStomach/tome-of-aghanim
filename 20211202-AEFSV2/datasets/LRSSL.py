import torch
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from datasets.base import Base

class LRSSL(Base):
    inited = False
    base = './data/LRSSL/'
    params = {
        'epoch': 1000, 'lr': 9e-05, 'wd': 6e-07,
        'sim_threshold': 0.5, 'loss_p_weight': 0.993, 'loss_d_weight': 0.994, 'loss_weight': 1,
        'a1': 0.00000000, 'a2': 0.00000000,
        
        'dropout': 0.2, 'graph_dropout': 0.1,
        'fr_dim': [4096, [1426, 1024], [682, 256], [4447, 2048]],
        'fd_dim': [2048, [1426, 1024], [682, 256], [4447, 2048]]
    }
    path = {
        # 'drugs': base + 'drug.txt',
        'drug_sim': base + 'drug_sim.txt',
        'drug_smiles': base + 'drug_smiles.csv',

        'drug_ecfps': base + 'drug_ecfps12.txt',
        'rpi': base + 'drug_target_domain_mat.txt', # 763*1426
        'rdi': base + 'drug_dis_mat.txt', # 763*682
    }

    def drugs(self):
        return np.loadtxt(self.path['rdi'], dtype=str, delimiter='\t')[1:, 0]

    def init(self, mask_drugs=None):
        self.mask_drugs = mask_drugs
        self.rpi = self.mask(self.data('rpi', skip=True))
        self.rdi = self.mask(self.data('rdi', skip=True))
        # self.rri = self.mask(self.data('rri'))

        drug_fps = self.mask(self.data('drug_ecfps', delimiter=','))
        # drug_vec = self.mask(self.data('drug_vec', delimiter=',')))
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter=',')
        ).T)

        self.drug_x1 = torch.from_numpy(drug_fps).float().to(self.device)
        self.drug_x2 = torch.from_numpy(np.matmul(self.drug_A, self.rpi)).float().to(self.device)
        self.drug_x3 = torch.from_numpy(np.matmul(self.drug_A, self.rdi)).float().to(self.device)
        # self.drug_x4 = torch.from_numpy(drug_vec).float().to(self.device)
        self.drug_edge = self.edge(self.drug_A, self.params['sim_threshold'])

        self.rnum = self.rpi.shape[0]
        self.pnum = self.rpi.shape[1]
        self.dnum = self.rdi.shape[1]

        self.inited = True
    
    def prepare(self):
        seqs = []
        radius = 6
        length = 4096

        drugs = np.loadtxt(self.path['drug_smiles'], delimiter=',', dtype=str, comments=None)
        for drug in drugs:
            try:
                name, smiles = drug
                mol = Chem.MolFromSmiles(smiles)
                seqs.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length).ToList())
            except Exception as e:
                print(drug, e)

        np.savetxt(self.path['drug_ecfps'], seqs, fmt='%s', delimiter=',')

    def data(self, name, dtype=int, delimiter='\t', skip=False):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return pd.read_csv(self.path[name], header=0, index_col=0, sep=delimiter).to_numpy(dtype) \
        if skip else \
        np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)
