import numpy as np
import json
import sys
import os
import math
import time
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
import joblib

os.system("taskset -p 0xff %d" % os.getpid())

threads = int(sys.argv[3])
#num_splits = 4
traindatafile = sys.argv[1]
outfile = sys.argv[2]

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


if __name__ == '__main__':
    print('Logistic regression')
    lr = LogisticRegression(C=1)
    lr.fit(X,y)
    if lr.classes_[0] == 1:
        poscol = 0
    else:
        poscol = 1
    joblib.dump(lr,outfile)
