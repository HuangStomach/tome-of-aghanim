import numpy as np
import pandas as pd
from grmf import GRMF
np.set_printoptions(suppress = True)
np.set_printoptions(threshold = np.inf)

mod = 'IC'
adj_mat = pd.read_table('./Data/adj_{}.txt'.format(mod), sep='\t', header=0, index_col=0).T
drug_mat = pd.read_table('./Data/drug_{}.txt'.format(mod), sep='\t', header=0, index_col=0)
target_mat = pd.read_table('./Data/target_{}.txt'.format(mod), sep='\t', header=0, index_col=0)

name = 'hsa55503'
index = 129
connect = adj_mat.loc[:, name].to_numpy().argsort()[::-1][:9]
adj_mat.loc[:, name] = adj_mat.loc[:, name].map(lambda x: 0)

grmf = GRMF()
grmf.fit(adj_mat, drug_mat, target_mat)
values = grmf.Y_hat[:, index]

print(connect)
print(np.argsort(values)[::-1][:20])
print(np.sort(values)[::-1][:20])
