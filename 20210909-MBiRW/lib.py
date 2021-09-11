import numpy as np
import pandas as pd
import random
import math
import os

def set_par_fun(ass_mat, sim_mat_s):
    '''
    ### Args:
        ass_mat: 药物和疾病的关联矩阵
        sim_mat: 药物或者疾病的相似矩阵
    '''
    sim_mat = sim_mat_s.copy()
    num = sim_mat.shape[0] # 样本数量
    s_mat = np.zeros(10)

    s = np.zeros(10, dtype=int)
    n = np.zeros(10, dtype=int)

    k = 0
    for i in range(1, num):
        for j in range(i):
            k += 1
            value = int(sim_mat[i][j] * 10)
            sharedValue = np.dot(ass_mat[i], ass_mat[j])
            if value > 9: value = 9

            s[value] += 1
            if sharedValue >= 1: n[value] += 1 # 样本和其他类别样本具有相似的作用样本

    for i in range(10):
        s_mat[i] = s[i] / n[i] # 依照样本相似度 转换成样本对儿的百分比
    
    random_mat = np.zeros((10, 10))
    # 使用洗牌法进行打乱
    for _i in range(10):
        s = np.zeros(10, dtype=int)
        n = np.zeros(10, dtype=int)

        m = 0
        rsum = int(num * (num - 1) / 2) # 所有sim_mat配对可能的数量
        dmat = np.zeros(rsum)

        # 打乱之后放回原处
        for i in range(1, num):
            for j in range(i):
                dmat[m] = sim_mat[i][j]
                m += 1
        
        for i in range(rsum):
            j = random.randrange(i, rsum)
            dmat[i], dmat[j] = dmat[j], dmat[i]
        
        m = 0
        for i in range(1, num):
            for j in range(i):
                sim_mat[i][j] = dmat[m];
                m += 1

        for i in range(1, num):
            for j in range(i):
                value = int(sim_mat[i][j] * 10)
                sharedValue = np.dot(ass_mat[i], ass_mat[j])
                if value > 9: value = 9

                s[value] += 1
                if sharedValue >= 1: 
                    n[value] += 1
        
        for i in range(10):
            random_mat[_i][i] = 0.0 if n[i] == 0 else s[i] / n[i]

    r_mat = np.mean(random_mat, axis=1) # 进行十次随机化结果，取平均值作为备用对照
    
    threshold = 0;
    for i in range(10):
        result_p = s_mat[i]
        result_p_random = r_mat[i]

        if result_p >= result_p_random:
            threshold = (i + 1) * 0.1;
            break

    return (math.log(99) - math.log(9999)) / threshold # 得到一个阈值来判断是否有参考价值

def cluster(drugs_sim, diseases_sim, drugs_share, diseases_share):
    '''
    ### Args:
        drugs_sim: 药物相似矩阵
        diseases_sim: 疾病相似矩阵
        drugs_share: 药物共享疾病个数矩阵
        diseases_share: 疾病共享药物个数矩阵
    '''
    drugs_name = pd.read_table('./Data/DrugsName', sep=' ', header=None, squeeze=True).to_numpy()
    diseases_name = pd.read_table('./Data/DiseasesName', sep=' ', header=None, squeeze=True).to_numpy()

    drugs_num = drugs_sim.shape[0]
    diseases_num = diseases_sim.shape[0]

    # 按照格式整理聚类模型所需文件
    file = open('./Output/drugsP.txt', 'w')
    for i in range(1, drugs_num):
        for j in range(i):
            shared_nums = drugs_share[i][j]
            if shared_nums > 0:
                content = "{}\t{}\t{:.0f}\n".format(drugs_name[i], drugs_name[j], shared_nums)
                file.write(content)
    file.close()
            
    file = open('./Output/diseasesP.txt', 'w')
    for i in range(1, diseases_num):
        for j in range(i):
            shared_nums = diseases_share[i][j]
            if shared_nums > 0:
                content = "{}\t{}\t{:.0f}\n".format(diseases_name[i], diseases_name[j], shared_nums)
                file.write(content)
    file.close()

    # 使用聚类模型进行聚类
    os.system('java -jar Lib/cluster_one-1.2.jar Output/drugsP.txt -F csv > Output/drugs_cluster.csv')
    os.system('java -jar Lib/cluster_one-1.2.jar Output/diseasesP.txt -F csv > Output/diseases_cluster.csv')

    # 记录
    drugs_cluster = pd.read_csv('./Output/drugs_cluster.csv')
    c_drugs = np.array([])
    c_drugs_quality = np.array([])
    for index, row in drugs_cluster.iterrows():
        pvalue = row['P-value']
        if pvalue > 0.1: continue
        
        quality = row['Quality']
        np.append(c_drugs, row['Members'])
        np.append(c_drugs_quality, quality)
    

    diseases_cluster = pd.read_csv('./Output/diseases_cluster.csv')
    c_diseases = np.array([])
    c_diseases_quality = np.array([])
    for index, row in diseases_cluster.iterrows():
        pvalue = row['P-value']
        if pvalue > 0.1: continue
        
        quality = row['Quality']
        np.append(c_diseases, row['Members'])
        np.append(c_diseases_quality, quality)

    drugs_cohesv = drugs_sim.copy()
    drugs_pos = np.zeros(100)
    for i in range(1, drugs_num):
        # 先查找所有药物，看每个药物是否存在于一个聚类中
        num = 0
        name = drugs_name[i]
        for j in range(c_drugs.shape[0]):
            if name not in c_drugs[j]: continue

            drugs_pos[num] = j
            num += 1
        
        # 如果只存在于1个聚类或根本没有聚类，则跳过
        if num <= 1: continue

        # 找到和其他药物内聚力最大的值作为参考记录
        for j in range(i):
            r_name = drugs_name[j]
            r_quality = 0
            flag = 0

            for k in range(num):
                index = drugs_pos[k]
                if r_name not in c_drugs[index]: continue

                flag = 1
                r_quality = c_drugs_quality[index] if r_quality < c_drugs_quality[index] else r_quality
            
            if flag != 1: continue
            drugs_cohesv[i][j] = (1 + r_quality) * drugs_sim[i][j]
            if drugs_cohesv[i][j] > 1: drugs_cohesv[i][j] = 0.99
            if drugs_cohesv[i][j] < r_quality: drugs_cohesv[i][j] = min(r_quality, 0.99)
            drugs_cohesv[j][i] = drugs_cohesv[i][j]
    for i in range(drugs_num): drugs_cohesv[i][i] = 1

    diseases_cohesv = diseases_sim.copy()
    diseases_pos = np.zeros(100)
    for i in range(1, diseases_num):
        # 先查找所有疾病，看每个疾病是否存在于一个聚类中
        num = 0
        name = diseases_name[i]
        for j in range(c_diseases.shape[0]):
            if name not in c_diseases[j]: continue

            diseases_pos[num] = j
            num += 1
        
        # 如果只存在于1个聚类或根本没有聚类，则跳过
        if num <= 1: continue

        # 找到和其他疾病内聚力最大的值作为参考记录
        for j in range(i):
            r_name = diseases_name[j]
            r_quality = 0
            flag = 0

            for k in range(num):
                index = diseases_pos[k]
                if r_name not in c_diseases[index]: continue

                flag = 1
                r_quality = c_diseases_quality[index] if r_quality < c_diseases_quality[index] else r_quality
            
            if flag != 1: continue
            diseases_cohesv[i][j] = (1 + r_quality) * diseases_sim[i][j]
            if diseases_cohesv[i][j] > 1: diseases_cohesv[i][j] = 0.99
            if diseases_cohesv[i][j] < r_quality: diseases_cohesv[i][j] = min(r_quality, 0.99)
            diseases_cohesv[j][i] = diseases_cohesv[i][j]
    for i in range(diseases_num): diseases_cohesv[i][i] = 1
            
    return (drugs_cohesv, diseases_cohesv)

def norm_fun(mat):
    num = mat.shape[0]
    sum_mat = np.zeros(num)
    result = np.zeros((num, num))

    for i in range(num):
        sum_mat[i] = np.sum(mat[i])

    for i in range(num):
        sum_a = sum_mat[i]
        for j in range(num):
            sum_b = sum_mat[j]
            if sum_a == 0 or sum_b == 0: result[i][j] = 0
            else: result[i][j] = mat[i][j] / math.sqrt(sum_a * sum_b)

    return result
