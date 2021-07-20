import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import joblib

# ccle模型pcc
models_info = joblib.load('../../Output/2/ccle_models_info.joblib')
tcga_pred_pd = pd.read_csv('../../Output/3/ccle_a_pred_tcga.csv')

plt.figure(figsize=(20, 8))
ax = plt.subplot()
ax.set_facecolor('whitesmoke')

ccle_mat = np.array([])
drugs_label = []
drugs_pos = []
i = 0
for drug, model_info in models_info.items():

    ccle_y = model_info['y'][model_info['y'] != -9]
    ccle_pred = model_info['pred'][model_info['pred'] != -9]
    tcga_pred = tcga_pred_pd.loc[:, drug]

    bplot = ax.boxplot([ccle_y, ccle_pred, tcga_pred], positions=[i + 0.3, i + 1, i + 1.7],
        widths=.7, patch_artist=True, sym='.k', medianprops=dict(color='k'))
    bplot['boxes'][0].set_facecolor('red')
    bplot['boxes'][1].set_facecolor('green')
    bplot['boxes'][2].set_facecolor('blue')

    i += 3
    drugs_pos.append(i)
    drugs_label.append(drug)

ax.set_xticks(range(1, len(drugs_label) * 3, 3))
ax.set_xticklabels(drugs_label, rotation = 45)

ax.plot([], c='red', label='CCLE Observed DR')
ax.plot([], c='green', label='CCLE Predicted DR')
ax.plot([], c='blue', label='TCGA Predicted DR')
ax.legend()
ax.set_ylabel('Drug response (ActArea)')

plt.show()
