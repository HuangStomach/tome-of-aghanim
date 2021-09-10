import numpy as np
import random
import math
import os

'''
ass_mat 药物和疾病的关联矩阵
sim_mat 药物或者疾病的相似矩阵
'''
def set_par_fun(ass_mat, sim_mat_s):
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

def cluster(drugs_sim, diseases_sim, drugs_share, diseases_share, drugs_name, diseases_name):
    drugs_num = drugs_sim.shape[0]
    diseases_num = diseases_sim.shape[0]

    file_drug = open('./Output/drugsP.txt', 'w')
    file_disease = open('./Output/diseasesP.txt', 'w')
    for i in range(1, drugs_num):
        for j in range(i):
            shared_nums = drugs_share[i][j]
            if shared_nums > 0:
                content = "{}\t{}\t{:.0f}\n".format(drugs_name[i], drugs_name[j], shared_nums)
                file_drug.write(content)
            
            shared_nums = diseases_share[i][j]
            if shared_nums > 0:
                content = "{}\t{}\t{:.0f}\n".format(diseases_name[i], diseases_name[j], shared_nums)
                file_disease.write(content)

    file_drug.close()
    file_disease.close()

    os.system('java -jar Lib/cluster_one-1.2.jar Output/drugsP.txt -F csv >> Output/drugs_cluster.csv')
    os.system('java -jar Lib/cluster_one-1.2.jar Output/diseasesP.txt -F csv >> Output/diseases_cluster.csv')
            
    return (1, 2)