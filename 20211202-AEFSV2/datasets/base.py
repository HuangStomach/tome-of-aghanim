from abc import ABCMeta
from abc import abstractmethod
import numpy as np

class Base:
    __metaclass__= ABCMeta

    @abstractmethod
    def drugs(self):
        pass

    def mask(self, mat):
        if self.mask_drugs is None: return mat
        mat = np.delete(mat, self.mask_drugs, axis=0)
        return mat
