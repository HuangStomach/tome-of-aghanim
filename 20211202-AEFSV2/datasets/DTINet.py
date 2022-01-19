import numpy as np
import torch
from sklearn.preprocessing import minmax_scale

from datasets.base import Base

class DTINet(Base):
    inited = False
    base = './data/DTINet/'
    path = {
        'drugs': base + 'drug.txt',
        'drug_sim': base + 'Similarity_Matrix_Drugs.txt',
        # 'drug_sim': 'base + 'Similarity_Matrix_Drugs_s.txt',

        'drug_ecfps': base + 'drug_ecfps12.txt',
        'drug_se': base + 'mat_drug_se.txt',
        'drug_vec': base + 'drug_vec.txt',
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
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter='    ')
        ).T)

        # self.drug_x1 = np.matmul(self.drug_A, drug_fps) # ecfps
        self.drug_x1 = torch.from_numpy(drug_fps).float().to(self.device)
        self.drug_x2 = torch.from_numpy(np.matmul(self.drug_A, self.rpi)).float().to(self.device)
        self.drug_x3 = torch.from_numpy(np.matmul(self.drug_A, self.rdi)).float().to(self.device)

        self.rnum = self.rpi.shape[0]
        self.pnum = self.rpi.shape[1]
        self.dnum = self.rdi.shape[1]

        self.inited = True

    def mask(self, mat):
        if self.mask_drugs is None: return mat
        mat = np.delete(mat, self.mask_drugs, axis=0)
        return mat

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)
