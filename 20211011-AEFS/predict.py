import torch
import numpy as np
from lib import *


if __name__ ==  '__main__':
    test_drug_fps = np.loadtxt('dataset/test_fps.txt')
    test_x = torch.from_numpy(test_drug_fps).float()

    model = torch.load("result/model.pkl")
    model.eval()
    preDTI, preRDA = model(test_x)
    np.savetxt('result/y_pre_DPI.txt', preDTI.detach().numpy(), fmt='%f')
