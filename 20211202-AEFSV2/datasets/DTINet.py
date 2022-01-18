from datasets.base import Base
import numpy as np

class DTINet(Base):
    prepared = False
    base = './data/DTINet/'
    path = {
        'drugs': base + 'drug.txt',
        'drug_sim': base + 'Similarity_Matrix_Drugs.txt',
        # 'drug_sim': './data/DTINet/Similarity_Matrix_Drugs_s.txt',

        'drug_ecfps': './data/DTINet/drug_ecfps8.txt',
        # 'rpi': './data/DTINet/mat_drug_protein.txt',
        'rpi': './data/DTINet/mat_drug_protein_s.txt',
        # 'rri': './data/DTINet/mat_drug_drug.txt',
        'rdi': './data/DTINet/mat_drug_disease.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')

    def prepare(self, mask_drugs=None):
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

        self.prepared = True

    def mask(self, mat):
        if self.mask_drugs is None: return mat
        mat = np.delete(mat, self.mask_drugs, axis=0)
        return mat

    def data(self, name, dtype=int, delimiter=' '):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)

    def edge(self, sim_mat):
        l = sim_mat.shape[0]
        
        edge_index = [[], []]
        edge_wight = []
        for i in range(l):
            for j in range(i + 1, l):
                if sim_mat[i][j] < 0.5: continue
                edge_index[0].append(i)
                edge_index[1].append(j)
                edge_wight.append(sim_mat[i][j])

        return (np.array(edge_index), np.array(edge_wight))
