from importlib import import_module
from feature import output

# classes = [GAT, GCN, GIN]
# for i in classes:
#     for j in classes:
#         pass


while True:
    datasets = ['DTINet', 'AEFS']
    for i, dataset in enumerate(datasets):
        print("[{}] {}".format(i, dataset))
    
    str_in = input("请选择数据源: ");
    if str_in.isdigit():
        index = int(str_in)
        if index < 0 or index >= len(datasets): continue
    else: continue
    dataset = getattr(import_module('dataset'), datasets[index])()
    break

while True:
    print("[0] 随机切分数据:")
    print("[1] 训练数据:")
    str_in = input("请选择需要进行的操作: ");
    if not str_in.isdigit(): continue

    index = int(str_in)
    if index < 0 or index >= 2: continue

    if index == 0:
        dataset.split_data()
        print('\033[32m完成\033[0m')
    elif index == 1:
        models = ['GAT', 'GCN', 'GIN']

        for model_drug in models:
            for model_target in models:
                output(model_drug, model_target, dataset)