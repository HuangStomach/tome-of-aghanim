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
    name = name.split('.')[0]
    if name[0] == 'X':
        name = name[1:]
    return name
rownames = np.array(list(map(normal_name, ccle_latent.index.values)))

coach = pd.read_csv('./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
drugs = pd.unique(coach['Compound'])
drugs.sort()
intersect = np.intersect1d(rownames, coach['CCLE Cell Line Name']) # 取细胞系交集
coach = coach.loc[coach['CCLE Cell Line Name'].isin(intersect)]

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
        # 选择不同药物下有交集的数据，使用np24的ccle数据做标注
        coach_drug = coach.loc[coach['Compound'] == drug]

        indexes = []
        for cell_name in coach_drug['CCLE Cell Line Name'].to_numpy(dtype=str):
            indexes.append(np.argwhere(rownames == cell_name)[0][0])
        
        x = ccle_latent.iloc[indexes].to_numpy(dtype=np.float)
        y = coach_drug['ActArea'].to_numpy(dtype=np.float)

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
            
            model.save('./Output/2/ccle_models/{}.h5'.format(drug))

joblib.dump(models_info, './Output/2/ccle_nn_models_info.joblib')  

