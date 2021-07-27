import joblib
import pandas as pd
import numpy as np
np.set_printoptions(suppress=True)
import scipy.stats as st

import tensorflow.keras.backend as K
from tensorflow.keras.layers import Dense, InputLayer
from tensorflow.keras import Sequential
from sklearn.metrics import r2_score

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
ccle_aliases = np.array(list(ccle_anno['Aliases']))

gdsc_in = gdsc_anno.loc[gdsc_anno['COSMIC_ID'].isin(ccle_cosmic_id)]
gdsc_not_in = gdsc_anno.loc[gdsc_anno['COSMIC_ID'].isin(ccle_cosmic_id) == False]

def r2(y_true, y_pred):
    a = K.square(y_pred - y_true)
    b = K.sum(a)
    c = K.mean(y_true)
    d = K.square(y_true - c)
    e = K.sum(d)
    f = 1 - b/e
    return f

end = 2

models_info = dict()
models = dict()

for step in range(1, end + 1):
    ccle_latent = pd.read_table('./Output/1/{}.CCLE_latent.tsv'.format(step), 
        index_col = 0, float_precision='high')
    # ccle_latent.index = pd.Series(rownames)

    for drug in drugs:
        # 有交集的数据
        coach_in = gdsc_in.loc[gdsc_in['DRUG_NAME'] == drug]
        gdsc_cosmic_id = list(coach_in['COSMIC_ID'])
        index1 = [ccle_cosmic_id.index(x) if x in ccle_cosmic_id else None for x in gdsc_cosmic_id]
        ccle_achs = ccle_anno.iloc[index1, 0]
        index1 = [rownames.index(x) if x in rownames else None for x in ccle_achs]

        coach_in = coach_in.iloc[[i for i, value in enumerate(index1) if not value is None]]
        index1 = [x for x in index1 if x is not None] # COSMIC_ID有交集的索引
        y = -coach_in.loc[:, "LN_IC50"].to_numpy()

        # 无交集的数据
        coach_not_in = gdsc_not_in.loc[gdsc_not_in['DRUG_NAME'] == '(5Z)-7-Oxozeaenol']
        def normal_cell(name):
            return name.replace('-', '')
        cellnames = np.array(list(map(normal_cell, coach_not_in['CELL_LINE_NAME'])))
        name_index = np.where(np.in1d(cellnames, ccle_aliases))[0]
        name_index = pd.unique(np.concatenate((name_index, np.where(np.in1d(coach_not_in['CELL_LINE_NAME'], ccle_aliases))[0])))

        coach_not_in = coach_not_in.iloc[name_index, ]
        gdsc_cellnames = coach_not_in['CELL_LINE_NAME'].to_numpy()
        cellnames = cellnames[name_index]

        index2 = np.array([])
        for k in range(len(cellnames)):
            if gdsc_cellnames[k] in ccle_aliases:
                indexes = np.argwhere(ccle_aliases == gdsc_cellnames[k])[0]
            else:
                indexes = np.argwhere(ccle_aliases == cellnames[k])[0]
            index2 = np.append(index2, indexes)

        index2 = [rownames.index(x) if x in rownames else None for x in ccle_anno.iloc[index2, 0]]
        coach_not_in = coach_not_in.iloc[[i for i, value in enumerate(index2) if not value is None]]
        gdsc_cellnames = coach_not_in['CELL_LINE_NAME'].to_numpy()

        index2 = np.array([])
        for k in range(len(gdsc_cellnames)):
            if gdsc_cellnames[k] in ccle_aliases:
                indexes = np.argwhere(ccle_aliases == gdsc_cellnames[k])
            else:
                indexes = np.argwhere(ccle_aliases == gdsc_cellnames[k].replace('-', ''))
            index2 = np.append(index2, indexes)

        index2 = [rownames.index(x) if x in rownames else None for x in ccle_anno.iloc[index2, 0]]
        y = np.append(y, -coach_not_in.loc[:, "LN_IC50"].to_numpy())
        indexes = np.append(index1, index2)
        x = ccle_latent.iloc[indexes]

        model = Sequential()
        model.add(InputLayer((x.shape[1], )))
        model.add(Dense(100, activation='relu'))
        model.add(Dense(1))
        model.compile(loss='mse', optimizer='adam', metrics=[r2])
        model.fit(x, y, epochs=100, validation_split=0.2, verbose=0)

        pred = model.predict(x, verbose=1)[:, 0]
        score = r2_score(y, pred)
        print(drug, score)

        if drug not in models_info.keys() or score > models_info[drug]['r2_score']:
            full_y = np.full(ccle_latent.shape[0], -9, dtype=float)
            full_pred = np.full(ccle_latent.shape[0], -9, dtype=float)
            full_y[indexes] = y
            full_pred[indexes] = pred

            models_info[drug] = {
                'y': full_y,
                'pred': full_pred,
                'step': step,
                'size': len(y),
                'pcc': st.pearsonr(y, pred),
                'r2_score': score,
            }
            model.save('./Output/2/gdsc_models/{}.h5'.format(drug))

joblib.dump(models_info, './Output/2/gdsc_nn_models_info.joblib')
