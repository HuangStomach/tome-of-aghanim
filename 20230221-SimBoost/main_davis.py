import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import igraph
from sklearn.decomposition import NMF

import xgboost
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RandomizedSearchCV

target_gene_names = pd.read_csv("data/target_gene_names.txt", header=None, index_col=0)
drug_pubchemIDs = pd.read_csv("data/drug_PubChem_CIDs.txt", header=None, index_col=0) 

f1, f2 = open("data/target_gene_names.txt"), open("data/drug_PubChem_CIDs.txt")
id1, id2 = [], []
for i in f1: id1.append(i.strip("\n"))
for i in f2: id2.append(int(i.strip("\n")))

sim_targets = pd.read_csv("data/target-target_similarities_WS_normalized.txt", 
    sep=" ", header=None, names=id1)
sim_targets.index = id1

# sim_drugs = pd.read_csv("data/drug_sis.csv", sep=",", header=None, names=id2)
sim_drugs = pd.read_csv("data/drug-drug_similarities_2D.txt", sep=" ", header=None, names=id2)
sim_drugs.index = id2

sim_drugs_g = pd.read_csv("data/drug_sis.csv", sep=",", header=None, names=id2)
sim_drugs_g.index = id2

bindings = pd.read_csv("data/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt", sep=" ", header=None,
    names=id1)
bindings.index = id2

target_gene_names.sort_index(inplace=True)
drug_pubchemIDs.sort_index(inplace=True)
bindings.sort_index(inplace=True)
sim_targets.sort_index(inplace=True)
sim_drugs.sort_index(inplace=True)

sim_targets = sim_targets / 100
transformed_bindings = -np.log10(bindings / (10 ** 9))

# plt.hist(transformed_bindings.stack())

l = []
for i in id1:
    for j in id2:
        k = transformed_bindings.loc[j, i]
        l2 = [j, i, k]
        l.append(l2)
drug_target_binding = pd.DataFrame(l, columns=['Drug', 'Target', 'Binding_Val'])

# drug_target_binding
train_data = None
val_data = None
test_data = None

train_data, x_remain = train_test_split(drug_target_binding, test_size=0.3)
val_data, test_data = train_test_split(x_remain, test_size=0.3)

for i in sim_targets:
    target_gene_names.loc[i, 't_avg-sim'] = np.mean(sim_targets.loc[i, :])

target_bindingvals = {}
for i in range(len(train_data)):
    x = train_data.iloc[i, :]  # drug=x[0],target=x[1],binding_val=[2]
    if x[1] not in target_bindingvals:
        target_bindingvals[x[1]] = [x[2]]
    else:
        target_bindingvals[x[1]].append(x[2])
for i in target_bindingvals:
    target_gene_names.loc[i, 't_avg-binding'] = np.mean(target_bindingvals[i])

for i in sim_drugs:
    drug_pubchemIDs.loc[i, 'd_avg-sim'] = np.mean(sim_drugs.loc[i, :])

drug_bindingvals = {}
for i in range(len(train_data)):
    x = train_data.iloc[i, :]  # drug=x[0],target=x[1],binding_val=[2]
    if x[0] not in drug_bindingvals:
        drug_bindingvals[x[0]] = [x[2]]
    else:
        drug_bindingvals[x[0]].append(x[2])
for i in drug_bindingvals:
    drug_pubchemIDs.loc[i, 'd_avg-binding'] = np.mean(drug_bindingvals[i])

sim_drugs = sim_drugs.reindex(sorted(sim_drugs.columns), axis=1)  # sorting columns
# sim_drugs = sim_drugs_g.reindex(sorted(sim_drugs_g.columns), axis=1)  # sorting columns

drug_sim_threshold = 0.6
target_sim_threshold = 0.6

drug_graph = igraph.Graph()
target_graph = igraph.Graph()

drug_graph.add_vertices(len(sim_drugs))
target_graph.add_vertices(len(sim_targets))

for i, drug_1 in enumerate(sim_drugs):
    for j, drug_2 in enumerate(sim_drugs):
        if (sim_drugs.loc[drug_1, drug_2] > drug_sim_threshold) and (drug_1 != drug_2):
            drug_graph.add_edges([(i, j)])

for i, tar_1 in enumerate(sim_targets):
    for j, tar_2 in enumerate(sim_targets):
        if (sim_targets.loc[tar_1, tar_2] > target_sim_threshold) and (tar_1 != tar_2):
            target_graph.add_edges([(i, j)])

for vertex in target_graph.vs:
    target_gene_names.loc[sorted(id1)[
        vertex.index], 't_n_neighbors'] = target_graph.neighborhood_size(vertex, mindist=1)
    target_gene_names.loc[sorted(id1)[vertex.index],
                          't_page_rank'] = target_graph.pagerank(vertex)


for vertex in drug_graph.vs:
    drug_pubchemIDs.loc[sorted(id2)[
        vertex.index], 'd_n_neighbors'] = drug_graph.neighborhood_size(vertex, mindist=1)
    drug_pubchemIDs.loc[sorted(id2)[vertex.index],
                        'd_page_rank'] = drug_graph.pagerank(vertex)

latent_dim = 3
train_binding_matrix = None
train_binding_matrix = pd.DataFrame(columns=sorted(id1), index=sorted(id2))  # empty dataframe

for i in range(len(train_data)):
    x = train_data.iloc[i, :]  # drug=x[0],target=x[1],binding_val=x[2]
    train_binding_matrix.loc[x[0], x[1]] = x[2]

train_binding_matrix = train_binding_matrix.fillna(5)

model = NMF(n_components=latent_dim, init='random', n_iter_int=1000, random_state=0)
P = model.fit_transform(train_binding_matrix)
Q = model.components_

for idx, x in enumerate(sorted(id2)):
    drug_pubchemIDs.loc[x, "d_features1"] = P[idx][0]
    drug_pubchemIDs.loc[x, "d_features2"] = P[idx][1]
    drug_pubchemIDs.loc[x, "d_features3"] = P[idx][2]

for idx, x in enumerate(sorted(id1)):
    target_gene_names.loc[x, "t_features1"] = Q.T[idx][0]
    target_gene_names.loc[x, "t_features2"] = Q.T[idx][1]
    target_gene_names.loc[x, "t_features3"] = Q.T[idx][2]

rows_list = []
for i in id1:  # id1=list of target names
    for j in id2:  # id2=list of drug names
        d_row = drug_pubchemIDs.loc[j, :].tolist()
        t_row = target_gene_names.loc[i, :].tolist()
        rows_list.append(t_row+d_row)

X = pd.DataFrame(rows_list, columns=[
    'd_avg-sim', 'd_avg-binding', 'd_n_neighbors', 
    'd_page_rank', 'd_features1', 'd_features2',
    'd_features3', 't_avg-sim', 't_avg-binding', 
    't_n_neighbors', 't_page_rank', 't_features1', 
    't_features2', 't_features3'
])
Y = drug_target_binding.loc[:, "Binding_Val"]

train_ratio = 0.7
val_ratio = 0.1
test_ratio = 0.2

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=1 - train_ratio)

X_val, X_test, Y_val, Y_test = train_test_split(
    X_test, Y_test, test_size=test_ratio/(test_ratio + val_ratio))

model = xgboost.XGBRegressor()

param_grid = {
    'learning_rate': np.arange(0.0005, 0.3, 0.0005),
    'n_estimators': [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 200],
    'max_depth': [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20],
    'colsample_bytree': np.arange(0.1, 1.0, 0.01),
    'subsample': np.arange(0.01, 1.0, 0.01)}

grid_search = RandomizedSearchCV(
    model, param_grid, random_state=0, cv=5, n_iter=10)
grid_result = grid_search.fit(X_train, Y_train)

def plot_model_results(results):
    epochs = len(results['validation_0']['rmse'])
    x_axis = range(0, epochs)
    fig, ax = plt.subplots()
    ax.plot(x_axis, results['validation_0']['rmse'], label='Train')
    ax.plot(x_axis, results['validation_1']['rmse'], label='Validation')
    ax.legend()
    plt.ylabel('RMSE')
    plt.show()


learning_rate = grid_result.best_params_['learning_rate']
n_estimators = grid_result.best_params_['n_estimators']
max_depth = grid_result.best_params_['max_depth']
colsample_bytree = grid_result.best_params_['colsample_bytree']
subsample = grid_result.best_params_['subsample']

model = xgboost.XGBRegressor(objective='reg:squarederror', learning_rate=learning_rate,
                             colsample_bytree=colsample_bytree,
                             max_depth=max_depth,
                             subsample=subsample,
                             n_estimators=n_estimators,
                             eval_metric='rmse')

model.fit(X_train, Y_train, eval_set=[(X_train, Y_train), (X_test, Y_test)], verbose=False)

validation_mse = mean_squared_error(Y_val, model.predict(X_val))
print("Validation MSE: %.3f" % validation_mse)
