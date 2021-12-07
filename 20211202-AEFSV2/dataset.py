import numpy as np
import sklearn
from lib import *

class Dataset:
    prepared = False
    path = {
        'drugs': './data/DTINet/drug.txt',
        'drug_sim': './data/DTINet/Similarity_Matrix_Drugs.txt',
        'protein_sim': './data/DTINet/Similarity_Matrix_Proteins.txt',

        'drug_fps': './data/DTINet/drug_ecfps.txt',
        'protein_embed': './data/DTINet/protein_embeds.csv',

        'rpi': './data/DTINet/mat_drug_protein.txt',
        'rri': './data/DTINet/mat_drug_drug.txt',
        'ppi': './data/DTINet/mat_protein_protein.txt',
        'rdi': './data/DTINet/mat_drug_disease.txt',
        'pdi': './data/DTINet/mat_protein_disease.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')

    def prepare(self, mask_drugs=None):
        print("Loading Data...")
        self.mask_drugs = mask_drugs
        self.rpi = self.mask(self.data('rpi'))
        self.rdi = self.mask(self.data('rdi'))
        self.pdi = self.data('pdi')
        self.ppi = self.data('ppi')

        drug_fps = self.mask(self.data('drug_fps', delimiter=','))
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter='    ')
        ).T)

        self.drug_x1 = np.matmul(self.drug_A, drug_fps) # ecfps
        self.drug_x2 = np.matmul(self.drug_A, self.rdi)
        self.drug_x3 = np.matmul(self.drug_A, self.rpi)

        self.drug_z1 = np.matmul(self.drug_A, self.rdi)
        self.drug_z2 = np.matmul(self.drug_A, self.rpi)

        self.rnum = drug_fps.shape[0]
        self.pnum = self.pdi.shape[0]
        self.dnum = self.pdi.shape[1]

        self.prepared = True

    def mask(self, mat):
        mat = mat[self.mask_drugs, :] = 0
        return mat

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)

    def edge(self, edge_mat, sim_mat):
        l = edge_mat.shape[0]
        
        edge_index = [[], []]
        edge_wight = []
        for i in range(l):
            for j in range(i + 1, l):
                if edge_mat[i][j] < 0.5: continue
                edge_index[0].append(i)
                edge_index[1].append(j)
                edge_wight.append(sim_mat[i][j])

        return (np.array(edge_index), np.array(edge_wight))

if __name__=='__main__':
    dataset = Dataset()
    dataset.prepare()