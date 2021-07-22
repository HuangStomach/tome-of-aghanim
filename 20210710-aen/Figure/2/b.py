import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

drugs_match_2 = pd.read_table("../../Data/Match/drugs_match_2.txt", header = None)
coach = pd.read_csv('../../Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
ccle_latent = pd.read_table('../../Output/1/1.CCLE_latent.tsv', index_col = 0)
def normal_name(name):
    name = name.split('.')[0]
    if name[0] == 'X':
        name = name[1:]
    return name
rownames = np.array(list(map(normal_name, ccle_latent.index.values)))

for match in drugs_match_2.to_numpy():
    coach_drug = coach.loc[coach['Compound'] == match[0]]
    shared_samples, shared_coach_index, shared_latent_index = np.intersect1d(coach_drug.iloc[:, 0].to_numpy(), rownames, return_indices=True)
    coach_drug = coach.iloc[shared_coach_index]
    ccle_latent_drug = ccle_latent.iloc[shared_latent_index]

    ccle_y = coach_drug.loc[:, "ActArea"]
    ccle_y_ic50 = coach_drug.loc[:, "IC50..uM."]