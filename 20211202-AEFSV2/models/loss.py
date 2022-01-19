import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable


class FocalLoss(nn.Module):
    def __init__(self, gamma=0, alpha=None, size_average=True):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.alpha = alpha
        if isinstance(alpha, (float, int)):
            self.alpha = torch.Tensor([alpha, 1-alpha])
        if isinstance(alpha, list):
            self.alpha = torch.Tensor(alpha)
        self.size_average = size_average

    def forward(self, input, target):
        if input.dim() > 2:
            # N,C,H,W => N,C,H*W
            input = input.view(input.size(0), input.size(1), -1)
            input = input.transpose(1, 2)    # N,C,H*W => N,H*W,C
            input = input.contiguous().view(-1, input.size(2))   # N,H*W,C => N*H*W,C
        target = target.view(-1, 1)

        logpt = F.log_softmax(input)
        logpt = input.gather(1, target)
        logpt = logpt.view(-1)
        pt = Variable(logpt.data.exp())

        if self.alpha is not None:
            if self.alpha.type() != input.data.type():
                self.alpha = self.alpha.type_as(input.data)
            at = self.alpha.gather(0, target.data.view(-1))
            logpt = logpt * Variable(at)

        loss = -1 * (1-pt)**self.gamma * logpt
        if self.size_average:
            return loss.mean()
        else:
            return loss.sum()

class WeightMSELoss(nn.Module):
    def __init__(self, alpha = 0.5):
        super(WeightMSELoss, self).__init__()
        self.alpha = alpha

    def forward(self, input, target):
        target_p = target
        target_n = torch.logical_not(target_p, out=torch.empty(target_p.size(), dtype=int, device=target_p.device))
        return self.alpha * F.mse_loss(target_p.mul(input), target) \
            + (1 - self.alpha) * F.mse_loss(target_n.mul(input), target)

class SONLoss(nn.Module):
    def __init__(self, k):
        '''
        ## Same order neighbours
        ### Parameters
            * k: top k neighbours
        '''
        super(SONLoss, self).__init__()
        self.k = k

    def forward(self, S_hat, S, eye):
        '''
        S: 药似性矩阵
        S: 计算出的相似性估计矩阵
        eye: 单位阵
        a: 约束参数
        '''
        _S_hat = S_hat - eye
        _S = S - eye
        S_hat_val, S_hat_idx = _S_hat.topk(self.k)
        S_val, S_idx = _S.topk(self.k)

        diff_idx = S_hat_idx.bitwise_xor(S_idx)
        cal_flag = diff_idx.div(diff_idx).nan_to_num(0)  # 索引相同为1 不同为0

        return S_hat_val.sub(S_val).pow(2).mul(cal_flag).sum(1).sqrt().sum()
