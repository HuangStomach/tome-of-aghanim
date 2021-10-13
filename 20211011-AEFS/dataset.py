import numpy as np

class Dataset:
    pass

class DTINet(Dataset):
    path = {
        'drugs': './datasets/DTINet/drug.txt',
        'drugs_fps': './datasets/DTINet/drug_ecfps.txt',
        'drugs_sim': './datasets/DTINet/Similarity_Matrix_Drugs.txt',
        'protein_sim': './datasets/DTINet/Similarity_Matrix_Proteins.txt',
        'dpi': './datasets/DTINet/mat_drug_protein.txt',
        'rda': './datasets/DTINet/mat_drug_disease.txt',
        'train_id': './datasets/DTINet/train_id.txt',
        'test_id': './datasets/DTINet/test_id.txt',
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
        return np.loadtxt(self.path['split_drug_sim'].format(type), dtype=float, delimiter='    ')
    
    def protein_sim(self):
        '''
        药物相似性
        '''
        return np.loadtxt(self.path['protein_sim'], dtype=float, delimiter='    ')
    
    def fps(self, type='train'):
        '''
        药物指纹
        '''
        return np.loadtxt(self.path['split_fps'].format(type), dtype=int, delimiter=' ')
    
    def dpi(self, type='train'):
        '''
        药物指纹
        '''
        return np.loadtxt(self.path['split_dpi'].format(type), dtype=int, delimiter=' ')
    
    def rda(self, type='train'):
        '''
        药物指纹
        '''
        return np.loadtxt(self.path['split_fps'].format(type), dtype=int, delimiter=' ')

    def split_data(self):
        drugs = np.loadtxt(self.path['drugs'], dtype=str, delimiter='\n')
        length = len(drugs)

        cut = int(length * 0.8)

        ids = np.arange(length)
        np.random.shuffle(ids)
        train_id, test_id = ids[:cut], ids[cut:]
        train_id.sort()
        test_id.sort()
        np.savetxt('./datasets/DTINet/train_id.txt', train_id, delimiter='\n', fmt='%d')
        np.savetxt('./datasets/DTINet/test_id.txt', test_id, delimiter='\n', fmt='%d')

        drug_sim = np.loadtxt(self.path['drugs_sim'], dtype=float, delimiter='    ')
        np.savetxt('./datasets/DTINet/train_drug_sim.txt', drug_sim[train_id], delimiter=' ', fmt='%d')
        np.savetxt('./datasets/DTINet/test_drug_sim.txt', drug_sim[test_id], delimiter=' ', fmt='%d')

        fps = np.loadtxt(self.path['drugs_fps'], dtype=int, delimiter=',')
        np.savetxt('./datasets/DTINet/train_fps.txt', fps[train_id], delimiter=' ', fmt='%d')
        np.savetxt('./datasets/DTINet/test_fps.txt', fps[test_id], delimiter=' ', fmt='%d')

        dpi = np.loadtxt(self.path['dpi'], dtype=int, delimiter=' ')
        np.savetxt('./datasets/DTINet/train_dpi.txt', dpi[train_id], delimiter=' ', fmt='%d')
        np.savetxt('./datasets/DTINet/test_dpi.txt', dpi[test_id], delimiter=' ', fmt='%d')

        rda = np.loadtxt(self.path['rda'], dtype=int, delimiter=' ')
        np.savetxt('./datasets/DTINet/train_rda.txt', rda[train_id], delimiter=' ', fmt='%d')
        np.savetxt('./datasets/DTINet/test_rda.txt', rda[test_id], delimiter=' ', fmt='%d')
