import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import joblib

# ccle模型pcc
models_info = joblib.load('../../Output/2/ccle_models_info.joblib')
ccle_pcc = np.array([])
for drug, model_info in models_info.items():
    info = np.array([drug, model_info['pcc'][0]])
    ccle_pcc = np.vstack((ccle_pcc, info)) if ccle_pcc.shape[0] > 1 else info

# gdsc模型pcc
models_info = joblib.load('../../Output/2/gdsc_models_info.joblib')
gdsc_pcc = np.array([])
for drug, model_info in models_info.items():
    info = np.array([drug, model_info['pcc'][0]])
    gdsc_pcc = np.vstack((gdsc_pcc, info)) if gdsc_pcc.shape[0] > 1 else info

drugs_match = pd.read_table("../../Data/Match/drugs_match_2.txt", header=None)
ccle_drugs = ccle_pcc[:, 0].tolist()
gdsc_drugs = gdsc_pcc[:, 0].tolist()
match_ccle = [ccle_drugs.index(x) if x in ccle_drugs else None for x in drugs_match.iloc[:, 0].to_list()]
match_gdsc = [gdsc_drugs.index(x) if x in gdsc_drugs else None for x in drugs_match.iloc[:, 2].to_list()]

plt.figure(figsize=(20, 10))
ax = plt.subplot()

x1 = np.delete(ccle_pcc, match_ccle, axis=0)
x1.sort(axis=0)
rect = Rectangle((0, 0), len(x1), 1, fc='#FFFFE0')
ax.add_patch(rect)
x1 = np.column_stack((x1, np.full(x1.shape[0], 0)))

x2 = np.delete(gdsc_pcc, match_gdsc, axis=0)
x2.sort(axis=0)
rect = Rectangle((len(x1), 0), len(x2), 1, fc='#E0FFFF')
ax.add_patch(rect)
x2 = np.column_stack((x2, np.full(x2.shape[0], 1)))
x2 = np.vstack((x1, x2))

x3 = ccle_pcc[match_ccle]
x3.sort(axis=0)
rect = Rectangle((len(x2), 0), len(x3), 1, fc='#90EE90')
ax.add_patch(rect)
x3 = np.column_stack((x3, np.full(x3.shape[0], 2)))
x3 = np.vstack((x2, x3))

x4 = gdsc_pcc[match_gdsc]
x4.sort(axis=0)
rect = Rectangle((len(x3), 0), len(x4), 1, fc='#ADD8E6')
ax.add_patch(rect)
x4 = np.column_stack((x4, np.full(x4.shape[0], 3)))
x4 = np.vstack((x3, x4))

x = list(range(len(x4)))
y = np.array(x4[:, 1], dtype=float)

ax.set_xlim([0, len(x4)])
ax.set_ylim([0, 1])
ax.set_yticks([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1])
ax.plot([0, 275], [0.5, 0.5])
ax.plot(x, y, '.', c='#555')
ax.set_ylabel('In-sample PCC')
ax.set_xlabel('Compound index')

for i, val in enumerate(y):
    if val >= 0.5 or x4[i][2] == '1': continue
    ax.text(i, val, x4[i][0])

plt.savefig('a.png')
