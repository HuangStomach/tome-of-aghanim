import numpy as np
from lib import *

class Dataset:
    def prepare(self):
        print("读取数据")

        self.drug_A = self.data('drug_sim', dtype=float, delimiter='    ')
        self.drug_x1 = np.matmul(self.drug_A, self.data('drug_fps', delimiter=',')) # 指纹
        self.drug_x2 = np.matmul(self.drug_A, self.data('rdi'))
        self.drug_x3 = np.matmul(self.drug_A, self.data('rpi'))

        self.protein_A = self.data('protein_sim', dtype=float)
        self.protein_x1 = self.data('protein_embed', dtype=float, delimiter=',') # word embedding
        self.protein_x2 = self.data('pdi')
        self.protein_x3 = self.data('rpi').T

    def edge_index(self, sim_mat):
        l = sim_mat.shape[0]
        
        edge_index = [[], []]
        for i in range(l):
            for j in range(i + 1, l):
                if sim_mat[i, j] < 0.5: continue
                edge_index[0].append(i)
                edge_index[1].append(j)
        return np.array(edge_index)

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, name):
            return getattr(self, name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)

class DTINet(Dataset):
    path = {
        'drugs': './datasets/DTINet/drug.txt',
        'drug_sim': './datasets/DTINet/Similarity_Matrix_Drugs.txt',
        'protein_sim': './datasets/DTINet/Similarity_Matrix_Proteins.txt',

        'drug_fps': './datasets/DTINet/drug_ecfps.txt',
        'protein_embed': './datasets/DTINet/protein_embeds.csv',

        'rpi': './datasets/DTINet/mat_drug_protein.txt',
        'rri': './datasets/DTINet/mat_drug_drug.txt',
        'ppi': './datasets/DTINet/mat_protein_protein.txt',
        'rdi': './datasets/DTINet/mat_drug_disease.txt',
        'pdi': './datasets/DTINet/mat_protein_disease.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
