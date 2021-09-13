import numpy as np
import pandas as pd
from lib import KNKN
np.set_printoptions(suppress = True)

mod = 'GPCR'
miu = 0.7
k = 5

adj_mat = pd.read_table('./Data/adj_{}.txt'.format(mod), sep='\t')
drug_mat = pd.read_table('./Data/drug_{}.txt'.format(mod), sep='\t')
target_mat = pd.read_table('./Data/target_{}.txt'.format(mod), sep='\t')
(m, n) = adj_mat.shape

Yd = drug_mat.copy()
Yt = target_mat.copy()
knkn = KNKN(k)
knkn.fit(drug_mat)
for d in range(n):
    (indexes, values) = knkn.neighbors(d)
    weights = np.zeros(k)
    for i in range(k):
        weights[i] = miu ** (i - 1) * values[i]
