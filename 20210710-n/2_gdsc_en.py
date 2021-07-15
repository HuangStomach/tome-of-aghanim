import joblib
import pandas as pd
import numpy as np
np.set_printoptions(suppress=True)
from sklearn.linear_model import ElasticNetCV
import scipy.stats as st

# ccle_to_tissue = np.load('./Output/1/ccle_to_tissue.npy', allow_pickle=True).item()
# tissues = np.unique(list(ccle_to_tissue.values()))
ccle_latent = pd.read_table('./Output/1/1.CCLE_latent.tsv', index_col = 0)

def normal_name(name):
    names = name.split('.')
    return '{}-{}'.format(names[2], names[3])
rownames = list(map(normal_name, ccle_latent.index.values))

ccle_anno = pd.read_csv("./Data/CCLE/DepMap-2018q3-celllines.csv", dtype={'COSMIC_ID': str})
gdsc_anno = pd.read_table("./Data/GDSC/v17.3_fitted_dose_response.txt", dtype={'COSMIC_ID': str})
drugs = pd.unique(gdsc_anno['DRUG_NAME'])
drugs.sort()

ccle_cosmic_id = list(ccle_anno['COSMIC_ID'])
gdsc_in = gdsc_anno.loc[gdsc_anno['COSMIC_ID'].isin(ccle_cosmic_id)]
gdsc_in = gdsc_in.loc[gdsc_anno['DRUG_NAME'] == '(5Z)-7-Oxozeaenol']
gdsc_not_in = gdsc_anno.loc[gdsc_anno['COSMIC_ID'].isin(ccle_cosmic_id) == False]

gdsc_cosmic_id = list(gdsc_in['COSMIC_ID'])
index1 = [ccle_cosmic_id.index(x) if x in ccle_cosmic_id else None for x in gdsc_cosmic_id]
ccle_achs = ccle_anno.iloc[index1, 0]
index1 = [rownames.index(x) if x in rownames else None for x in ccle_achs]
index1 = list(filter(None, index1))
gdsc_in = gdsc_in.iloc[index1]

end = 1

models_info = dict()
models = dict()

for step in range(1, end + 1):
    ccle_latent = pd.read_table('./Output/1/{}.CCLE_latent.tsv'.format(step), 
        index_col = 0, float_precision='high')
    # ccle_latent.index = pd.Series(rownames)

    for drug in drugs:
        # 选择不同药物下有交集的数据，使用np24的ccle数据做标注
        coach_drug = gdsc_in.loc[gdsc_anno['DRUG_NAME'] == drug]
        y = -coach_drug.loc[:, "LN_IC50"]

        indexes = []
        for cell_name in coach_drug['CCLE Cell Line Name'].to_numpy(dtype=str):
            indexes.append(np.argwhere(rownames == cell_name)[0][0])
        
        ccle_latent_drug = ccle_latent.iloc[indexes].to_numpy(dtype=np.float)
        y = coach_drug['ActArea'].to_numpy(dtype=np.float)

        regr = ElasticNetCV(cv=10, max_iter=10000, random_state=0)
        regr.fit(ccle_latent_drug, y)
        pred = regr.predict(ccle_latent_drug)
        score = regr.score(ccle_latent_drug, y)

        if drug not in models_info.keys() or score > models_info[drug]['r2_score']:
            models_info[drug] = {
                'size': len(y),
                'pcc': st.pearsonr(y, pred),
                'r2_score': score,
            }
            models[drug] = regr
            continue

joblib.dump(models_info, './Output/2/ccle_models_info.joblib')  
joblib.dump(models, './Output/2/ccle_models.joblib')  
