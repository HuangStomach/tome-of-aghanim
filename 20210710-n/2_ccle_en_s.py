import pandas as pd
import numpy as np

coach = pd.read_csv('./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
drugs = pd.unique(coach['Compound'])
drugs.sort()

ccle_latent = pd.read_table('./Output/1/1.CCLE_latent.tsv', index_col = 0)

def normal_name(name):
    name = name.split('.')[0]
    if name[0] == 'X':
        name = name[1:]
    return name
rownames = np.array(list(map(normal_name, ccle_latent.index.values)))
cell_names = rownames[np.char.find(rownames, 'HAEMATOPOIETIC') == -1] # 对细胞系名称进行格式化 去除造血细胞
intersect = np.intersect1d(cell_names, coach['CCLE Cell Line Name']) # 取细胞系交集

coach = coach.loc[coach['CCLE Cell Line Name'].isin(intersect)]

end = 1
for step in range(1, end + 1):
    ccle_latent = pd.read_table('./Output/1/{}.CCLE_latent.tsv'.format(step), index_col = 0)
    ccle_latent.index = pd.Series(rownames)

    for drug in drugs:
        coach_drug = coach.loc[coach['Compound'] == drug]
        indexs = [i for i, name in enumerate(rownames) if name in coach_drug['CCLE Cell Line Name'].to_numpy()]
        ccle_latent_drug = ccle_latent.iloc[indexs]
        print(ccle_latent_drug)