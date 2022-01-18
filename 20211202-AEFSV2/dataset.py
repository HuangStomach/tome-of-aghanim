import numpy as np
import importlib

class Dataset:
    def __init__(self, type="DTINet"):
        module = importlib.import_module('datasets.{}'.format(type))
        self.handler = getattr(module, type)()

    def __getattr__(self, name):
        return getattr(self.handler, name)

    def drugs(self):
        return self.handler.drugs()

    def prepare(self, mask_drugs=None):
        print("Loading Data...")
        return self.handler.prepare(mask_drugs)

    def mask(self, mat):
        return self.handler.mask(mat)

    def data(self, name, dtype=int, delimiter=' '):
        return self.handler.data(name, dtype, delimiter)

    def edge(self, sim_mat, threshold):
        return self.handler.edge(sim_mat, threshold)

if __name__=='__main__':
    dataset = Dataset()
    dataset.prepare()
