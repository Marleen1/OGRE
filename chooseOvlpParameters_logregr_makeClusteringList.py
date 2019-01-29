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
from sklearn import svm
from sklearn.neural_network import MLPClassifier as MLP
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier as RFC
import multiprocessing as mp
from scipy import stats

os.system("taskset -p 0xff %d" % os.getpid())

threads = int(sys.argv[4])
#num_splits = 4
ovlpfile = sys.argv[1]
traindatafile = sys.argv[2]
outfile = sys.argv[3]

'''
data_same = []
with open(datafile_same,'r') as f:
    for line in f:
        data_same.append(line.rstrip().split('\t'))
data_same = np.vstack([np.asarray(row)[2:6].astype(float) for row in data_same])

data_s = np.hstack((data_same,np.ones((data_same.shape[0],1)))).astype(float)
data_same = None
print(data_s.shape)

data_diff = []
with open(datafile_diff,'r') as f:
    for line in f:
        data_diff.append(line.rstrip().split('\t'))
data_diff = np.vstack([np.asarray(row)[2:6].astype(float) for row in data_diff])

data_d = np.hstack((data_diff,np.zeros((data_diff.shape[0],1)))).astype(float)
data_diff = None
print(data_d.shape)

data = np.vstack((data_s,data_d))
'''

data_str = np.genfromtxt(traindatafile,delimiter='\t',dtype=str)
data = data_str[:,2:].astype(float)

#===============
# Learning
#===============

ncols = data.shape[1]
X = data[:,:ncols-1]
y = data[:,ncols-1]
ncols = X.shape[1]

Xnew = np.hstack((np.expand_dims(X[:,0],axis=1),np.expand_dims(X[:,2],axis=1)))
ncols = Xnew.shape[1]
X = Xnew
#X = Xnew-np.amin(Xnew)
#X = X/np.amax(X)

'''
Xnew = X
for i in range(ncols):
    for j in range(i+1,ncols):
        newcol = np.multiply(X[:,i],X[:,j])
        Xnew = np.hstack((Xnew,np.expand_dims(newcol,axis=1)))
X = Xnew
'''

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
    print('Logistic regression')
    lr = LogisticRegression(C=1)
    lr.fit(X,y)
    if lr.classes_[0] == 1:
        poscol = 0
    else:
        poscol = 1
    print('Split files')
    num_lines = sum(1 for line in open(ovlpfile))
    lines_per_file = 10000000 #math.ceil(num_lines/num_splits)
    num_splits = math.ceil(num_lines/10000000)
    if num_splits > 900:
        lines_per_file = math.ceil(num_lines/900)
        num_splits = 900
    print(num_splits)
    os.system('split -l '+str(lines_per_file)+' -a 3 --numeric-suffixes=1 --additional-suffix=.txt '+ovlpfile+' ovlp_split')
    print('Process files')
    writetofile('001')
    po = mp.Pool(threads)
    results = po.map(writetofile, [str(p).zfill(3) for p in range(1,num_splits+1)])
    po.close()
    po.join()
#    accs = np.empty((num_splits,1),dtype=float)
#    tps = np.empty((num_splits,1),dtype=float)
#    tns = np.empty((num_splits,1),dtype=float)
#    fps = np.empty((num_splits,1),dtype=float)
#    fns = np.empty((num_splits,1),dtype=float)
#    for i in range(num_splits):
#        accs[i] = results[i][0]
#        tns[i] = results[i][1]
#        fps[i] = results[i][2]
#        fns[i] = results[i][3]
#        tps[i] = results[i][4]
#    with open('accuracy_'+outfile,'w') as f:
#        f.write('Validation accuracy: '+str(np.mean(accs))+'\n')
#        f.write('tps: '+str(np.sum(tps))+'\n')
#        f.write('fps: '+str(np.sum(fps))+'\n')
#        f.write('tns: '+str(np.sum(tns))+'\n')
#        f.write('fns: '+str(np.sum(fns))+'\n')
    print('Concatenate files')
    os.system('cat ovlp_predict_keep_split* > ovlp_predict_temp.txt')
    os.system('cat ovlp_predict_removed_split* > ovlp_predict_removed.txt')
    os.system('/bin/rm ovlp_predict_keep_split*')
    os.system('/bin/rm ovlp_predict_removed_split*')
#    os.system('/bin/rm ovlp_split*')
    print('Sort')
    os.system('sort -nr -k3,3 ovlp_predict_temp.txt > '+outfile)
    os.system('/bin/rm ovlp_predict_temp.txt')
