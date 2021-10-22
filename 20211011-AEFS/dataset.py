import numpy as np
import torch
from lib import *

class Dataset:
    def prepare(self):
        print("读取数据")
        drug_A = self.drug_sim('train')
        drug_x1 = np.matmul(drug_A, self.fps('train')) # 指纹
        drug_x2 = np.matmul(drug_A, self.dpi('train'))
        drug_x3 = np.matmul(drug_A, self.rda('train'))

        protein_A = max_min_normalize(self.protein_sim())
        protein_x1 = np.matmul(protein_A, self.protein_embed()) # wordembedding
        protein_x2 = np.matmul(protein_A, self.dpi('train').T)
        protein_x3 = np.matmul(protein_A, self.pdi())
        print(drug_x1.shape, drug_x2.shape, drug_x3.shape)
        print(protein_x1.shape, protein_x2.shape, protein_x3.shape)
        
        print("numpy 转 tensor")
        self.id = torch.from_numpy(self.idx('train')).int()
        self.drug_x1 = torch.from_numpy(drug_x1).float()
        self.drug_x2 = torch.from_numpy(drug_x2).float()
        self.drug_x3 = torch.from_numpy(drug_x3).float()

        self.protein_x1 = torch.from_numpy(protein_x1).float()
        self.protein_x2 = torch.from_numpy(protein_x2).float()
        self.protein_x3 = torch.from_numpy(protein_x3).float()

        self.drug_A = torch.from_numpy(drug_A).float()
        self.protein_A = torch.from_numpy(protein_A).float()

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
        return np.loadtxt(self.path['split_fps'].format(type), dtype=int, delimiter=' ')
    
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
