import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fig, (ax0, ax1) = plt.subplots(1, 2)
fig.set_size_inches(14, 10)

tcga_pred_pd = pd.read_csv('../../Output/3/ccle_nn_pred_tcga.csv').loc[:, ['Sample', 'Cancer', 'Lapatinib']]
drugs = tcga_pred_pd.columns[2:]
cancer_types = pd.unique(tcga_pred_pd['Cancer'])
sample_type = np.array([x[13:15] for x in tcga_pred_pd.iloc[:, 0].to_numpy()], dtype=str)
# 取14，15位 14-15位为2位数字，01-09表示肿瘤样本，10-16表示正常对照样本

cancer_drug_ccle = np.array([])
for cancer in cancer_types:
    type_code = "01"
    if cancer == "LAML": type_code = "03"
    if cancer == "SKCM": type_code = "06"

    tmp_pd = tcga_pred_pd.iloc[np.where(sample_type == type_code)]
    tmp_pd = tmp_pd.loc[tmp_pd['Cancer'] == cancer]

    cancer_drug_ccle = np.vstack((cancer_drug_ccle, tmp_pd)) if cancer_drug_ccle.shape[0] > 1 else tmp_pd

clin_data = pd.read_table('./BRCA_clinicalMatrix')
clin_type = np.array([x[13:15] for x in clin_data.iloc[:, 0].to_numpy()], dtype=str)
clin_data = clin_data.iloc[np.where(clin_type == '01')]

status_name = 'lab_proc_her2_neu_immunohistochemistry_receptor_status' # 太长了 缩缩
state = clin_data.loc[clin_data[status_name].isin(["Equivocal", "Negative", "Positive"]), ["sampleID", status_name]]
clin_data = clin_data.iloc

samples_id = cancer_drug_ccle[:, 0].tolist()
x = state.iloc[np.where(state.iloc[:, 0].isin(samples_id))].iloc[:, 1]
indexes = [samples_id.index(x) if x in samples_id else None for x in state.iloc[:, 0]]
indexes = list(filter(None, indexes)) 
brca_ccle = cancer_drug_ccle[indexes]

data = pd.DataFrame(data={'Lapatinib': brca_ccle[:, 2], 'x': x})

ax0.set_facecolor('gainsboro')

dataset = [
    data[data['x'] == 'Equivocal']['Lapatinib'],
    data[data['x'] == 'Negative']['Lapatinib'],
    data[data['x'] == 'Positive']['Lapatinib']
]
bplot = ax0.boxplot(dataset,
    widths=.7, patch_artist=True, sym='', 
    medianprops=dict(color='k'), zorder=1,
    labels=["Equivocal", "Negative", "Positive"]
)

colors = ['red', 'green', 'blue']
for i, box in enumerate(bplot['boxes']):
    box.set_facecolor(colors[i])
    y = dataset[i]
    x = np.random.normal(i+1, .1, size=len(y))
    ax0.plot(x, y, '.k', markersize=2, zorder=2)
    ax0.text(i + .95, 1.5, len(y))

ax0.set_title('CCLE-based Model\np = ')
ax0.set_ylabel('Predicted Response to Lapatinib')

tcga_pred_pd = pd.read_csv('../../Output/3/gdsc_nn_pred_tcga.csv').loc[:, ['Sample', 'Cancer', 'Lapatinib']]
drugs = tcga_pred_pd.columns[2:]
cancer_types = pd.unique(tcga_pred_pd['Cancer'])
sample_type = np.array([x[13:15] for x in tcga_pred_pd.iloc[:, 0].to_numpy()], dtype=str)
# 取14，15位 14-15位为2位数字，01-09表示肿瘤样本，10-16表示正常对照样本

cancer_drug_gdsc = np.array([])
for cancer in cancer_types:
    type_code = "01"
    if cancer == "LAML": type_code = "03"
    if cancer == "SKCM": type_code = "06"

    tmp_pd = tcga_pred_pd.iloc[np.where(sample_type == type_code)]
    tmp_pd = tmp_pd.loc[tmp_pd['Cancer'] == cancer]

    cancer_drug_gdsc = np.vstack((cancer_drug_gdsc, tmp_pd)) if cancer_drug_gdsc.shape[0] > 1 else tmp_pd

clin_data = pd.read_table('./BRCA_clinicalMatrix')
clin_type = np.array([x[13:15] for x in clin_data.iloc[:, 0].to_numpy()], dtype=str)
clin_data = clin_data.iloc[np.where(clin_type == '01')]

status_name = 'lab_proc_her2_neu_immunohistochemistry_receptor_status' # 太长了 缩缩
state = clin_data.loc[clin_data[status_name].isin(["Equivocal", "Negative", "Positive"]), ["sampleID", status_name]]
clin_data = clin_data.iloc

samples_id = cancer_drug_gdsc[:, 0].tolist()
x = state.iloc[np.where(state.iloc[:, 0].isin(samples_id))].iloc[:, 1]
indexes = [samples_id.index(x) if x in samples_id else None for x in state.iloc[:, 0]]
indexes = list(filter(None, indexes)) 
brca_gdsc = cancer_drug_gdsc[indexes]

data = pd.DataFrame(data={'Lapatinib': brca_gdsc[:, 2], 'x': x})

ax1.set_facecolor('gainsboro')

dataset = [
    data[data['x'] == 'Equivocal']['Lapatinib'],
    data[data['x'] == 'Negative']['Lapatinib'],
    data[data['x'] == 'Positive']['Lapatinib']
]
bplot = ax1.boxplot(dataset,
    widths=.7, patch_artist=True, sym='', 
    medianprops=dict(color='k'), zorder=1,
    labels=["Equivocal", "Negative", "Positive"]
)

colors = ['red', 'green', 'blue']
for i, box in enumerate(bplot['boxes']):
    box.set_facecolor(colors[i])
    y = dataset[i]
    x = np.random.normal(i+1, .1, size=len(y))
    ax1.plot(x, y, '.k', markersize=2, zorder=2)
    ax1.text(i + .95, -1.3, len(y))

ax1.set_title('GDSC-based Model\np = ')
ax1.set_ylabel('Predicted Response to Lapatinib')

plt.savefig('a_nn.png', bbox_inches='tight')
