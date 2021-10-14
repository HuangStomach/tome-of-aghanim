import torch
from torch.nn import nn

# same order neighbours
class SONLoss(nn.Module):
    def __init__(self, Sr, k):
        super(SONLoss, self).__init__()
        self.Sr = Sr
        self.k = k
    
    def forward(self, S, a):
        '''
        Sr: 药物相似性矩阵
        S: 蛋白或疾病分数矩阵
        k: 邻居数
        a: 约束参数
        '''
        S_sum = S.pow(2).sum(dim=1)
        S1 = S_sum.unsqueeze(1)
        S2 = S_sum.unsqueeze(0)
        S_dist = torch.sqrt(S1 + S2 - 2 * S.mm(S.T)) # 与其他样本的距离

        S_val, S_idx = S_dist.topk(self.k + 1, largest=False)
        Sr_val, Sr_idx, Sr_val = self.Sr.topk(self.k + 1)
        S_val, S_idx, Sr_val, Sr_idx = S_val[:, 1:], S_idx[:, 1:], Sr_val[:, 1:], Sr_idx[:, 1:]

        diff_idx = S_idx.bitwise_xor(Sr_idx)
        punish = diff_idx.abs().sum(1) # 对位近邻的索引之差视为惩罚项
        cal_flag = diff_idx.div(diff_idx).nan_to_num(0) # 索引相同为1 不同为0
        
        return a * S_val.sub(Sr_val).pow(2).mul(cal_flag).sum(1).sqrt().mul(punish).sum()