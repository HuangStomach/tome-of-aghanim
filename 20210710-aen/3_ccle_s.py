import pandas as pd
import numpy as np
import joblib
models = joblib.load('./Output/2/ccle_s_models.joblib')
models_info = joblib.load('./Output/2/ccle_s_models_info.joblib')
tcga_to_cancer = joblib.load('./Output/1/tcga_to_cancer.joblib') 

# ccle模型预测tcga
tcga_pred = np.array([])
columnnames = np.array(['Cancer', 'Sample'])
for drug, model in models.items():
    tcga_latent = pd.read_table('./Output/1/{}.TCGA_latent.tsv'.format(models_info[drug]['step']), 
        index_col = 0, float_precision='high')
    result = model.predict(tcga_latent)
    
    tcga_pred = np.vstack((tcga_pred, result)) if tcga_pred.shape[0] > 1 else result
    columnnames = np.append(columnnames, drug)
    print('tcga:', drug)

def normal_name(name):
    return name.replace('.', '-')
rownames = np.array(list(map(normal_name, tcga_latent.index.values)))

def normal_cancer(name):
    return tcga_to_cancer[name]
cancers = np.array(list(map(normal_cancer, tcga_latent.index.values)))

tcga_pred = np.vstack((cancers, tcga_pred))
tcga_pred = np.vstack((rownames, tcga_pred))

tcga_pred_df = pd.DataFrame(tcga_pred.T, columns=columnnames)
tcga_pred_df.to_csv('./Output/3/ccle_s_pred_tcga.csv', index=False)

# ccle模型自我预测
ccle_latent = pd.read_table('./Output/1/1.CCLE_latent.tsv', index_col = 0)
rownames = ccle_latent.index.values
ccle_pred = np.array([])
columnnames = np.array(['Cellline'])
for drug, model_info in models_info.items():
    ccle_pred = np.vstack((ccle_pred, model_info['pred'])) if ccle_pred.shape[0] > 1 else model_info['pred']

ccle_pred = np.vstack((rownames, ccle_pred))
ccle_pred_df = pd.DataFrame(ccle_pred.T, columns=columnnames)
ccle_pred_df.to_csv('./Output/3/ccle_s_pred_ccle.csv', index=False)
