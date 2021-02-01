#!/usr/bin/env python

import pandas as pd
import json
import sys
import copy
import seaborn as sns
import matplotlib.pyplot as plt

clusters_json=sys.argv[1]
dist_file=sys.argv[2]
mash_size=sys.argv[3] #17000
out_file=sys.argv[4]


c2read=json.load(open(clusters_json,'r'))

c2c={} # original cluster id to merged cluster id, '1':set(1,2)  '2':set(1,2)
i=0
with open(dist_file,'r') as fr:
    for line in fr:
        i+=1
        #if i%1000==0:print("### "+str(i))
        ci,cj=line.split("\t")[0:2]


        seti=c2c[ci] if ci in c2c.keys() else set([ci])
        setj=c2c[cj] if cj in c2c.keys() else set([cj])
        reads=[]
        for s in seti|setj:
            reads.extend(c2read[s])

        if mash_size == 'inf':
            pass
        elif len(reads)>int(mash_size):
            continue

        for s in seti|setj:
            c2c[s]=seti|setj

cc=set([':'.join(list(val)) for key,val in c2c.items()])
merged_c={}
merged_c2read={}
for clusters in cc:
    clusters=clusters.split(':')
    new_key=clusters[0]
    new_val=[]
    for c in clusters:
        new_val.extend(c2read[c])
        merged_c[c]=1
    merged_c2read[new_key]=new_val

for k,v in c2read.items():
    if k not in merged_c.keys():
        merged_c2read[k]=v
with open(out_file,"w") as fw:
    json.dump(merged_c2read,fw)

