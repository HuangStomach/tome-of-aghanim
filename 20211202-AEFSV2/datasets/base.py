import torch
import numpy as np
from abc import ABCMeta
from abc import abstractmethod

class Base:
    __metaclass__= ABCMeta

    def __init__(self) -> None:
        self.device = 'cuda:0' if torch.cuda.is_available() else 'cpu'


    @abstractmethod
    def drugs(self):
        pass

    def mask(self, mat):
        if self.mask_drugs is None: return mat
        mat = np.delete(mat, self.mask_drugs, axis=0)
        return mat
    
    def edge(self, sim_mat, threshold = 0.5):
        l = sim_mat.shape[0]
        
        edge_index = [[], []]
        edge_wight = []
        for i in range(l):
            for j in range(i + 1, l):
                if sim_mat[i][j] < threshold: continue
                edge_index[0].append(i)
                edge_index[1].append(j)
                edge_wight.append(sim_mat[i][j])

        drug_edge = torch.tensor(edge_index).long().to(self.device)
        drug_weight = torch.tensor(edge_wight).long().to(self.device)
        return (drug_edge, drug_weight)
