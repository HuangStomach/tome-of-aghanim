import numpy as np
from lib import *

class Dataset:
    def prepare(self):
        print("读取数据")
        self.id = self.idx('train')
        self.drug_A = self.drug_sim('train')
        self.drug_x1 = np.matmul(self.drug_A, self.fps('train')) # 指纹
        self.drug_x2 = np.matmul(self.drug_A, self.rda('train'))
        self.drug_x3 = np.matmul(self.drug_A, self.dpi('train'))

        self.protein_x1 = self.protein_embed() # wordembedding
        self.protein_x2 = self.pdi()
        self.protein_x3 = self.dpi('train').T
    
    def edge_index(self, sim_mat):
        l = sim_mat.shape[0]
        edge_index = []
        for i in range(l):
            for j in range(i + 1, l):
                if sim_mat[i, j] >= 0.5:
                    edge_index.append([i, j])
        return np.array(edge_index)

class DTINet(Dataset):
    path = {
        'drugs': './datasets/DTINet/drug.txt',
        'drugs_fps': './datasets/DTINet/drug_ecfps.txt',
        'drugs_sim': './datasets/DTINet/Similarity_Matrix_Drugs.txt',
        'protein_sim': './datasets/DTINet/Similarity_Matrix_Proteins.txt',
        'protein_embed': './datasets/DTINet/protein_embeds.csv',
        'dpi': './datasets/DTINet/mat_drug_protein.txt',
        'rda': './datasets/DTINet/mat_drug_disease.txt',
        'pdi': './datasets/DTINet/mat_protein_disease.txt',
        'split_id': './datasets/DTINet/{}_id.txt',
        'split_drug_sim': './datasets/DTINet/{}_drug_sim.txt',
        'split_fps': './datasets/DTINet/{}_fps.txt',
        'split_dpi': './datasets/DTINet/{}_dpi.txt',
        'split_rda': './datasets/DTINet/{}_rda.txt',
    }

    def drugs(self):
        return np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
    
    def drug_sim(self, type='train'):
        '''
        药物相似性
        '''
        return np.loadtxt(self.path['split_drug_sim'].format(type), dtype=float, delimiter=' ')
    
    def idx(self, type='train'):
        '''
        药物index
        '''
        return np.loadtxt(self.path['split_id'].format(type), dtype=float, delimiter='\n')
    
    def protein_sim(self):
        '''
        药物相似性
        '''
        return np.loadtxt(self.path['protein_sim'], dtype=float, delimiter=' ')
    
    def protein_embed(self):
        '''
        药物相似性
        '''
        return np.loadtxt(self.path['protein_embed'], dtype=float, delimiter=',')
    
    def fps(self, type='train'):
        '''
        药物指纹
        '''
        return np.loadtxt(self.path['split_fps'].format(type), dtype=int, delimiter=' ')
    
    def dpi(self, type='train'):
        '''
        药物蛋白关联
        '''
        return np.loadtxt(self.path['split_dpi'].format(type), dtype=int, delimiter=' ')
    
    def rda(self, type='train'):
        '''
        药物疾病关联
        '''
        return np.loadtxt(self.path['split_rda'].format(type), dtype=int, delimiter=' ')
    
    def pdi(self):
        '''
        蛋白疾病关联
        '''
        return np.loadtxt(self.path['pdi'].format(type), dtype=int, delimiter=' ')

    def split_data(self):
        drugs = np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
        length = len(drugs)

        cut = int(length * 0.8)

        ids = np.arange(length)
        np.random.shuffle(ids)
        train_id, test_id = ids[:cut], ids[cut:]
        train_id.sort()
        test_id.sort()
        np.savetxt(self.path['split_id'].format('train'), train_id, delimiter='\n', fmt='%d')
        np.savetxt(self.path['split_id'].format('test'), test_id, delimiter='\n', fmt='%d')

        drug_sim = np.loadtxt(self.path['drugs_sim'], dtype=float, delimiter='    ')
        np.savetxt(self.path['split_drug_sim'].format('train'), drug_sim[train_id][:, train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_drug_sim'].format('test'), drug_sim[test_id][:, test_id], delimiter=' ', fmt='%d')

        fps = np.loadtxt(self.path['drugs_fps'], dtype=int, delimiter=',')
        np.savetxt(self.path['split_fps'].format('train'), fps[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_fps'].format('test'), fps[test_id], delimiter=' ', fmt='%d')

        dpi = np.loadtxt(self.path['dpi'], dtype=int, delimiter=' ')
        np.savetxt(self.path['split_dpi'].format('train'), dpi[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_dpi'].format('test'), dpi[test_id], delimiter=' ', fmt='%d')

        rda = np.loadtxt(self.path['rda'], dtype=int, delimiter=' ')
        np.savetxt(self.path['split_rda'].format('train'), rda[train_id], delimiter=' ', fmt='%d')
        np.savetxt(self.path['split_rda'].format('test'), rda[test_id], delimiter=' ', fmt='%d')

    # def split_graph(self):
    #     drug_sim = np.loadtxt(self.path['drugs_sim'], dtype=float, delimiter='    ')
    #     d = drug_sim[690]
    #     d.sort()
    #     print(d)
    #     quit()
    #     drugs = self.drugs()

    #     drug_mark = np.zeros(drugs.shape[0])

    #     graphs = []
    #     print(drugs.shape[0])
    #     for i in range(drugs.shape[0]):
    #         if drug_mark[i] != 0: continue

    #         drug_mark[i] = 1
    #         g = [i]
    #         self.dfs(drug_sim, drug_mark, i, g)
    #         graphs.append(g)
        
    #     for g in graphs:
    #         print(len(g))

    # def dfs(self, drug_sim, drug_mark, i, g):
    #     for j, drug in enumerate(drug_sim[i]):
    #         if drug >= 0.5 and drug_mark[j] == 0:
    #             g.append(j)
    #             drug_mark[j] = 1
    #             self.dfs(drug_sim, drug_mark, j, g)