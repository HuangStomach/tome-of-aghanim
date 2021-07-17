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
    name = name.split('.')[0]
    if name[0] == 'X':
        name = name[1:]
    return name
rownames = np.array(list(map(normal_name, ccle_latent.index.values)))
cell_names = rownames[np.char.find(rownames, 'HAEMATOPOIETIC') == -1] # 对细胞系名称进行格式化 去除造血细胞

coach = pd.read_csv('./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
drugs = pd.unique(coach['Compound'])
drugs.sort()
intersect = np.intersect1d(cell_names, coach['CCLE Cell Line Name']) # 取细胞系交集
coach = coach.loc[coach['CCLE Cell Line Name'].isin(intersect)]

end = 1

models_info = dict()
models = dict()

for step in range(1, end + 1):
    ccle_latent = pd.read_table('./Output/1/{}.CCLE_latent.tsv'.format(step), 
        index_col = 0, float_precision='high')
    # ccle_latent.index = pd.Series(rownames)

    for drug in drugs:
        # 选择不同药物下有交集的数据，使用np24的ccle数据做标注
        coach_drug = coach.loc[coach['Compound'] == drug]

        indexes = []
        for cell_name in coach_drug['CCLE Cell Line Name'].to_numpy(dtype=str):
            indexes.append(np.argwhere(rownames == cell_name)[0][0])
        
        ccle_latent_drug = ccle_latent.iloc[indexes].to_numpy(dtype=np.float)
        y = coach_drug['ActArea'].to_numpy(dtype=np.float)

        regr = ElasticNetCV(cv=10, max_iter=100000, random_state=0)
        regr.fit(ccle_latent_drug, y)
        pred = regr.predict(ccle_latent_drug)
        score = regr.score(ccle_latent_drug, y)
        print(drug, score)
        
        if drug not in models_info.keys() or score > models_info[drug]['r2_score']:
            models_info[drug] = {
                'y': y,
                'pred': pred,
                'step': step,
                'size': len(y),
                'pcc': st.pearsonr(y, pred),
                'r2_score': score,
            }
            models[drug] = regr
            continue

joblib.dump(models_info, './Output/2/ccle_s_models_info.joblib')  
joblib.dump(models, './Output/2/ccle_s_models.joblib')  
