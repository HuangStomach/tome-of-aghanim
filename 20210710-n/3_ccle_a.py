import pandas as pd
import numpy as np
import joblib
models = joblib.load('./Output/2/ccle_models.joblib')
models_info = joblib.load('./Output/2/ccle_models_info.joblib')
tcga_to_cancer = joblib.load('./Output/1/tcga_to_cancer.joblib') 

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
tcga_pred_df.to_csv('./Output/3/ccle_a_pred_tcga.csv', index=False)
