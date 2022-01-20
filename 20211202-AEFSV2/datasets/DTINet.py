import numpy as np
import torch
from sklearn.preprocessing import minmax_scale

from datasets.base import Base

class DTINet(Base):
    inited = False
    base = './data/DTINet/'
    params = {
        'epoch': 1000, 'lr': 9e-05, 'wd': 6e-07,

        'sim_threshold': 0.5, 'loss_p_weight': 0.998, 'loss_d_weight': 0.95, 'loss_weight': 0.0001,

        'a1': 0.000000001, 'a2': 0.000000001,
    }
    path = {
        'drugs': base + 'drug.txt',
        'drug_sim': base + 'Similarity_Matrix_Drugs.txt',

        'drug_ecfps': base + 'drug_ecfps12.txt',
        'rpi': base + 'mat_drug_protein.txt',
        # 'rpi': base + 'mat_drug_protein_s.txt',
        'rri': base + 'mat_drug_drug.txt',
        'rdi': base + 'mat_drug_disease.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')

    def init(self, mask_drugs=None):
        self.mask_drugs = mask_drugs
        self.rpi = self.mask(self.data('rpi'))
        self.rdi = self.mask(self.data('rdi'))
        # self.rri = self.mask(self.data('rri'))

        drug_fps = self.mask(self.data('drug_ecfps', delimiter=','))
        # drug_vec = self.mask(self.data('drug_vec', delimiter=',')))
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter='    ')
        ).T)

        # self.drug_x1 = np.matmul(self.drug_A, drug_fps) # ecfps
        self.drug_x1 = torch.from_numpy(drug_fps).float().to(self.device)
        self.drug_x2 = torch.from_numpy(np.matmul(self.drug_A, self.rpi)).float().to(self.device)
        self.drug_x3 = torch.from_numpy(np.matmul(self.drug_A, self.rdi)).float().to(self.device)
        # self.drug_x4 = torch.from_numpy(drug_vec).float().to(self.device)

        self.rnum = self.rpi.shape[0]
        self.pnum = self.rpi.shape[1]
        self.dnum = self.rdi.shape[1]

        self.inited = True

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)
