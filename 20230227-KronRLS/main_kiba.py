# http://staff.cs.utu.fi/~aatapa/software/RLScore/tutorial_kernels.html#tutorial-3-kronecker-learners
import numpy as np
from rlscore.learner import KronRLS
from sklearn.metrics import mean_squared_error

def load_davis(sis):
    Y = np.loadtxt("./data/kiba/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt")
    Y = -np.log10(Y / (10 ** 9))
    if sis:
        XD = np.loadtxt("./data/kiba/drug_sis.csv", delimiter=",")
    else:
        XD = np.loadtxt("./data/kiba/drug-drug_similarities_2D.txt")
    XT = np.loadtxt("./data/kiba/target-target_similarities_WS_normalized.txt")    
    return XD, XT, Y

def settingD_split(sis=False):
    np.random.seed(1)
    XD, XT, Y = load_davis(sis)
    drug_ind = list(range(Y.shape[0]))
    target_ind = list(range(Y.shape[1]))
    np.random.shuffle(drug_ind)
    np.random.shuffle(target_ind)
    train_drug_ind = drug_ind[:40]
    test_drug_ind = drug_ind[40:]
    train_target_ind = target_ind[:300]
    test_target_ind = target_ind[300:]
    #Setting 4: ensure that d,t pairs do not overlap between
    #training and test set
    Y_train = Y[np.ix_(train_drug_ind, train_target_ind)]
    Y_test = Y[np.ix_(test_drug_ind, test_target_ind)]
    Y_train = Y_train.ravel(order='F')
    Y_test = Y_test.ravel(order='F')
    XD_train = XD[train_drug_ind]
    XT_train = XT[train_target_ind]
    XD_test = XD[test_drug_ind]
    XT_test = XT[test_target_ind]
    return XD_train, XT_train, Y_train, XD_test, XT_test, Y_test

if __name__=="__main__":
    X1_train, X2_train, Y_train, X1_test, X2_test, Y_test = settingD_split(False)
    learner = KronRLS(X1 = X1_train, X2 = X2_train, Y = Y_train)
    log_regparams = np.arange(15, 35)
    for log_regparam in log_regparams:
        learner.solve(2.**log_regparam)
        P = learner.predict(X1_test, X2_test)
        mse = mean_squared_error(Y_test, P)
        print("regparam 2**%d, mse %f" %(log_regparam, mse))
