import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from lib import *

class LRSSL:
    inited = False
    base = './data/LRSSL/'
    path = {
        # 'drugs': base + 'drug.txt',
        'drug_sim': base + 'Similarity_Matrix_Drugs.txt',
        'drug_smiles': base + 'drug_smiles.csv',

        'drug_ecfps': base + 'drug_ecfps12.txt',
        'rpi': base + 'mat_drug_protein.txt',
        'rri': base + 'mat_drug_drug.txt',
        'rdi': base + 'drug_dis_mat.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['rdi'], dtype=str, delimiter='\t')[1:, 0]

    def init(self, mask_drugs=None):
        print("Loading Data...")
        self.mask_drugs = mask_drugs
        self.rpi = self.mask(self.data('rpi'))
        self.rdi = self.mask(self.data('rdi'))

        drug_fps = self.mask(self.data('drug_ecfps', delimiter=','))
        # drug_se = self.mask(self.data('drug_se'))
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter='    ')
        ).T)

        # self.drug_x1 = np.matmul(self.drug_A, drug_fps) # ecfps
        self.drug_x1 = drug_fps
        self.drug_x2 = np.matmul(self.drug_A, self.rdi)
        self.drug_x3 = np.matmul(self.drug_A, self.rpi)
        # self.drug_x4 = drug_se

        self.drug_z1 = np.matmul(self.drug_A, self.rdi)
        self.drug_z2 = np.matmul(self.drug_A, self.rpi)

        self.rnum = drug_fps.shape[0]
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
        
    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)
