import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

tcga_pred_pd = pd.read_csv('../../Output/3/ccle_nn_pred_tcga.csv').loc[:, ['Sample', 'Cancer', 'Paclitaxel']]
cancer_types = pd.unique(tcga_pred_pd['Cancer'])
sample_type = np.array([x[13:15] for x in tcga_pred_pd.iloc[:, 0].to_numpy()], dtype=str)

cancer_drug_ccle = np.array([])
for cancer in cancer_types:
    type_code = "01"
    if cancer == "LAML": type_code = "03"
    if cancer == "SKCM": type_code = "06"

    tmp_pd = tcga_pred_pd.iloc[np.where(sample_type == type_code)]
    tmp_pd = tmp_pd.loc[tmp_pd['Cancer'] == cancer]

    cancer_drug_ccle = np.vstack((cancer_drug_ccle, tmp_pd)) if cancer_drug_ccle.shape[0] > 1 else tmp_pd

sample_name = [x[0:12] for x in cancer_drug_ccle[:, 0]]
response = pd.read_table('../../Data/Response/drug_response.txt')
response = response.loc[response['drug.name'] == 'Paclitaxel'].to_numpy()
indexes = [sample_name.index(x) if x in sample_name else None for x in response[:, 1]]
indexes = np.array([[i, x] for i, x in enumerate(indexes) if x is not None])

response = np.hstack((response[indexes[:, 0]], cancer_drug_ccle[indexes[:, 1]]))
#[cancers patient.arr drug.name drug.id response start.time end.time pathological_time procurement_method procurement_time nte_time nte_pharm nte_surgery, sample_name, cancer, Paclitaxel]
indexes = [i for i, x in enumerate(response) if 'Disease' in x[4]]
p = stats.ttest_ind(response[indexes, -1], np.delete(response[:, -1], indexes)).pvalue

plt.figure(figsize=(8, 10))
ax = plt.subplot()
ax.set_facecolor('whitesmoke')
labels = ["Clinical Progressive Disease", "Stable Disease", "Partial Response", 'Complete Response']

dataset = [
    response[response[:, 4] == labels[0]][:, -1],
    response[response[:, 4] == labels[1]][:, -1],
    response[response[:, 4] == labels[2]][:, -1],
    response[response[:, 4] == labels[3]][:, -1],
]
bplot = ax.boxplot(dataset,
    widths=.7, patch_artist=True, sym='', 
    medianprops=dict(color='k'), zorder=1,
    labels=labels
)

colors = ['lightcoral', 'yellowgreen', 'mediumturquoise', 'plum']
for i, box in enumerate(bplot['boxes']):
    box.set_facecolor(colors[i])
    y = dataset[i]
    x = np.random.normal(i + 1, .1, size=len(y))
    ax.plot(x, y, '.k', markersize=2, zorder=2)
    ax.text(i + .95, 7.25, len(y))

ax.set_title('CCLE-based Model\np = {}'.format(p))
ax.set_ylabel('Predicted Response to Paclitaxel, TCGA-BRCA')

plt.savefig('d_nn.png', bbox_inches='tight')
