import numpy as np
from sklearn import metrics
import importlib

class Dataset:
    def __init__(self, type="DTINet"):
        module = importlib.import_module('datasets.{}'.format(type))
        self.handler = getattr(module, type)()

    def __getattr__(self, name):
        return getattr(self.handler, name)

    def drugs(self):
        return self.handler.drugs()

    def init(self, mask_drugs=None):
        print("Loading Data...")
        return self.handler.init(mask_drugs)

    def prepare(self):
        return self.handler.prepare()

    def data(self, name, dtype=int, delimiter=' '):
        return self.handler.data(name, dtype, delimiter)

    def splits(self):
        if hasattr(self.handler, 'splits'): self.handler.splits()
        
        drug_count = self.drugs().shape[0]
        shuffled_drugs = np.arange(drug_count)
        np.random.shuffle(shuffled_drugs)
        return np.array_split(shuffled_drugs, 10)

    def metric(self, target, target_hat):
        aupr_list = []
        for i, row in enumerate(target):
            if np.sum(row) == 0 : continue
            aupr_list.append(metrics.average_precision_score(row, target_hat[i]))

        target_f = target.flatten()
        target_hat_f = target_hat.flatten()
        fpr, tpr, _ = metrics.roc_curve(target_f, target_hat_f, pos_label=1)
        auc = metrics.auc(fpr, tpr)
        aupr = metrics.average_precision_score(target_f, target_hat_f)
        aupr_mean = np.mean(aupr_list)

        return (auc, aupr, aupr_mean)

if __name__=='__main__':
    dataset = Dataset()
    dataset.init()
