import numpy as np
import json
import sys
import os
import math
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from scipy import stats
import joblib

os.system("taskset -p 0xff %d" % os.getpid())

regrmodel = sys.argv[1]
splitnum = sys.argv[2]


def writetofile(fileID):
    data = np.genfromtxt('ovlp_split'+fileID+'.txt',delimiter='\t',dtype=str)
    readIDs = data[:,0:2]
    X = data[:,2:6].astype(float)
    Xnew = np.hstack((np.expand_dims(X[:,0],axis=1),np.expand_dims(X[:,2],axis=1)))
    ncols = Xnew.shape[1]
    X = Xnew
#    y_true = data[:,-1].astype(int)
#    acc = lr.score(X,y_true)
    y = lr.predict_proba(X)
#    y_hat = lr.predict(X)
#    tn, fp, fn, tp = confusion_matrix(y_true,y_hat).ravel()
    data_out_all = np.hstack((readIDs,np.expand_dims(y[:,poscol],axis=1)))
    data_out = data_out_all[np.where(data_out_all[:,2].astype(float) >= 0.5)]
    data_removed = data_out_all[np.where(data_out_all[:,2].astype(float) < 0.5)]
    data_out_all = None
    np.savetxt('ovlp_predict_keep_split'+fileID+'.txt',data_out,delimiter='\t',fmt='%s')
    np.savetxt('ovlp_predict_removed_split'+fileID+'.txt',data_removed,delimiter='\t',fmt='%s')
    return
#    return [acc, tn, fp, fn, tp]

if __name__ == '__main__':
    lr = joblib.load(regrmodel)
    if lr.classes_[0] == 1:
        poscol = 0
    else:
        poscol = 1
    writetofile(splitnum)
