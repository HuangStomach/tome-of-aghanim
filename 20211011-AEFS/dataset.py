import numpy as np
from lib import *

class Dataset:
    def prepare(self):
        print("读取数据")

        self.id = self.split('id', delimiter='\n')
        self.drug_A = self.split('drug_sim', dtype=float)
        self.drug_x1 = np.matmul(self.drug_A, self.split('fps')) # 指纹
        self.drug_x2 = np.matmul(self.drug_A, self.split('rdi'))
        self.drug_x3 = np.matmul(self.drug_A, self.split('rpi'))

        self.protein_A = self.protein_sim()
        self.protein_x1 = self.protein_embed() # word embedding
        self.protein_x2 = self.split('pdi')
        self.protein_x3 = self.split('rpi').T

    def edge_index(self, sim_mat):
        l = sim_mat.shape[0]
        
        edge_index = [[], []]
        for i in range(l):
            for j in range(i + 1, l):
                if sim_mat[i, j] < 0.5: continue
                edge_index[0].append(i)
                edge_index[1].append(j)
        return np.array(edge_index)

    def split(self, name, type='train', dtype=int, delimiter=' '):
        if hasattr(self, name):
            return getattr(self, name)()

        key = "split_{}".format(name)
        if key not in self.path: return []

        return np.loadtxt(self.path[key].format(type), dtype=dtype, delimiter=delimiter)

class DTINet(Dataset):
    path = {
        'drugs': './datasets/DTINet/drug.txt',
        'drugs_fps': './datasets/DTINet/drug_ecfps.txt',
        'drugs_sim': './datasets/DTINet/Similarity_Matrix_Drugs.txt',
        'protein_sim': './datasets/DTINet/Similarity_Matrix_Proteins.txt',
        'protein_embed': './datasets/DTINet/protein_embeds.csv',

        'rpi': './datasets/DTINet/mat_drug_protein.txt',
        'rri': './datasets/DTINet/mat_drug_drug.txt',
        'ppi': './datasets/DTINet/mat_protein_protein.txt',
        'rdi': './datasets/DTINet/mat_drug_disease.txt',
        'pdi': './datasets/DTINet/mat_protein_disease.txt',

        'split_id': './datasets/DTINet/{}_id.txt',
        'split_drug_sim': './datasets/DTINet/{}_drug_sim.txt',
        'split_rri': './datasets/DTINet/{}_rri.txt',
        'split_fps': './datasets/DTINet/{}_fps.txt',
        'split_rpi': './datasets/DTINet/{}_rpi.txt',
        'split_rdi': './datasets/DTINet/{}_rdi.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
    
    def protein_sim(self):
        '''
        蛋白相似性
        '''
        return np.loadtxt(self.path['protein_sim'], dtype=float, delimiter=' ')
    
    def ppi(self):
        '''
        蛋白的关系
        '''
        return np.loadtxt(self.path['ppi'], dtype=int, delimiter=' ')

    def pdi(self):
        '''
        蛋白疾病关联
        '''
        return np.loadtxt(self.path['pdi'], dtype=int, delimiter=' ')
    
    def protein_embed(self):
        '''
        蛋白质embedding
        '''
        return np.loadtxt(self.path['protein_embed'], dtype=float, delimiter=',')

    def split_data(self):
        drugs = np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
        length = len(drugs)

        cut = int(length * 0.05)

        ids = np.arange(length)
        np.random.shuffle(ids)
        train_id, test_id = ids[:cut], ids[cut:]
        train_id.sort()
        test_id.sort()
        np.savetxt(self.path['split_id'].format('train'), train_id, delimiter='\n', fmt='%d')
        np.savetxt(self.path['split_id'].format('test'), test_id, delimiter='\n', fmt='%d')

        drug_sim = np.loadtxt(self.path['drugs_sim'], dtype=float, delimiter='    ')
        np.savetxt(self.path['split_drug_sim'].format('train'), drug_sim[train_id][:, train_id], delimiter=' ')
        np.savetxt(self.path['split_drug_sim'].format('test'), drug_sim[test_id][:, test_id], delimiter=' ')

        rri = np.loadtxt(self.path['rri'], dtype=int, delimiter=' ')
        np.savetxt(self.path['split_rri'].format('train'), rri[train_id][:, train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_rri'].format('test'), rri[test_id][:, test_id], delimiter=' ', fmt='%d')

        fps = np.loadtxt(self.path['drugs_fps'], dtype=int, delimiter=',')
        np.savetxt(self.path['split_fps'].format('train'), fps[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_fps'].format('test'), fps[test_id], delimiter=' ', fmt='%d')

        rpi = np.loadtxt(self.path['rpi'], dtype=int, delimiter=' ')
        np.savetxt(self.path['split_rpi'].format('train'), rpi[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_rpi'].format('test'), rpi[test_id], delimiter=' ', fmt='%d')

        rdi = np.loadtxt(self.path['rdi'], dtype=int, delimiter=' ')
        np.savetxt(self.path['split_rdi'].format('train'), rdi[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_rdi'].format('test'), rdi[test_id], delimiter=' ', fmt='%d')
