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

    def split(self):
        rpi = self.data('rpi')
        rdis = []
        for i in range(rpi.shape[0]):
            for j in range(rpi.shape[1]):
                if rpi[i, j] == 1:
                    rdis.append((i, j))
        return rdis

    def prepare(self, ignore_drugs=None):
        print("Loading Data...")
        self.rpi = self.data('rpi')
        self.rdi = self.data('rdi')
        self.pdi = self.data('pdi')
        self.rri = self.data('rri')
        self.ppi = self.data('ppi')

        drug_fps = self.data('drug_fps', delimiter=',')
        self.drug_A = self.data('drug_sim', dtype=float, delimiter='    ')
        
        protein_embed = self.data('protein_embed', dtype=float, delimiter=',')
        self.protein_A = sklearn.preprocessing.minmax_scale(self.data('protein_sim', dtype=float))

        self.drug_x1 = np.matmul(self.drug_A, drug_fps) # ecfps
        self.drug_x2 = np.matmul(self.drug_A, self.rdi)
        self.drug_x3 = np.matmul(self.drug_A, self.rpi)

        self.protein_x1 = np.matmul(self.protein_A, protein_embed) # lstm embedding
        self.protein_x2 = np.matmul(self.protein_A, self.pdi)

        self.rnum = drug_fps.shape[0]
        self.pnum = protein_embed.shape[0]
        self.dnum = self.pdi.shape[1]

        self.prepared = True

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

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)

if __name__=='__main__':
    dataset = Dataset()
    dataset.prepare()