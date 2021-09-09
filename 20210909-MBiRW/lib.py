import numpy as np
import random
import math

'''
ass_mat 药物和疾病的关联矩阵
sim_mat 药物或者疾病的相似矩阵
'''
def set_par_fun(ass_mat, sim_mat):
    num = sim_mat.shape[0]
    s_mat = np.zeros(10)

    s = np.zeros(10)
    n = np.zeros(10)

    for i in range(1, num):
        for j in range(i + 1):
            value = int(sim_mat[i][j] * 10)
            sharedValue = np.dot(ass_mat[i], ass_mat[j])
            if value > 9: value = 9

            s[value] += 1
            if sharedValue >= 1: n[value] += 1

    for i in range(10):
        s_mat[i] = s[i] / n[i]
    
    result_mat = np.zeros((10, 10))
    for _i in range(10):
        # 使用洗牌法进行打乱
        s = np.zeros(10)
        n = np.zeros(10)

        m = 0
        rsum = num * (num - 1) / 2 # 所有sim_mat配对可能的数量
        dmat = np.zeros((1, rsum))

        # 打乱之后放回原处
        for i in range(1, num):
            for j in range(i + 1):
                dmat[m] = sim_mat[i][j]
                m += 1
        
        for i in range(rsum):
            j = random.randrange(i, rsum)
            dmat[i], dmat[j] = dmat[j], dmat[i]
        
        m = 0
        for i in range(1, num):
            for j in range(i + 1):
                sim_mat[i][j] = dmat[m];
                m += 1

        for i in range(1, num):
            for j in range(i + 1):
                value = int(sim_mat[i][j] * 10)
                sharedValue = np.dot(ass_mat[i], ass_mat[j])
                if value > 9: value = 9

                s[value] += 1
                if sharedValue >= 1: n[value] += 1

    for i in range(10):
        result_mat[_i][i] = s[i] / n[i]
    r_mat = np.mean(result_mat, axis=1)
    
    threshold = 0;
    for i in range(10):
        result_p = s_mat[i]
        result_pr = r_mat[i]

        if result_p >= 1 * result_pr:
            threshold = i * 0.1;
            break

    return (math.log(99) - math.log(9999)) / threshold
