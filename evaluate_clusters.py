import json
import numpy as np
import operator as op
import sys
from functools import reduce
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import math

dataset = sys.argv[1]
data_clustersizes = sys.argv[2]
datafolder = sys.argv[3]

#------------
# Statistics
#------------

clusterlist1 = json.load(open(data_clustersizes,'r'))
clusterlist = dict((k,v) for k,v in clusterlist1.items() if int(v) >= 20)
clusters = json.load(open(dataset+'_clusters_grouped.json','r')) #'/linuxhome/tmp/marleen/Minimap/cluster/limited/until30/minimap_clusterdict_until30_limited.json','r'))

'''
i = 0
for cl, reads in clusters.items():
    print(cl)
    print(reads)
    for rd in [reads]:
        print(rd)
    i += 1
    if i > 10:
        break

for k,v in clusterlist.items():
    print(k)
    print(v)
    break

sys.exit()

i = 0
for k,v in clusters.items():
    print(k)
    print(v)
    i += 1
    if i > 10:
        break
'''

with open('cluster_evaluation_file_'+dataset+'.txt','w') as f:
    f.write('Number of clusters: '+str(len(clusterlist.keys()))+'\n')
    f.write('Min cluster size of clusters: '+str(np.amin([val for val in clusterlist.values()]))+'\n')
    f.write('Mean cluster size of clusters: '+str(np.mean([val for val in clusterlist.values()]))+'\n')
    f.write('Max cluster size of clusters: '+str(np.amax([val for val in clusterlist.values()]))+'\n')

    f.write('Number of clusters with single read: '+str(len([k for k, v in clusterlist.items() if v==1]))+'\n')
    f.write('10th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],10))+'\n')
    f.write('20th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],20))+'\n')
    f.write('30th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],30))+'\n')
    f.write('40th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],40))+'\n')
    f.write('50th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],50))+'\n')
    f.write('60th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],60))+'\n')
    f.write('70th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],70))+'\n')
    f.write('80th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],80))+'\n')
    f.write('90th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],90))+'\n')
    f.write('92nd percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],92))+'\n')
    f.write('95th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],95))+'\n')
    f.write('97th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],97))+'\n')
    f.write('99th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],99))+'\n')
    f.write('100th percentile of clustersizes: '+str(np.percentile([v for v in clusterlist.values() if v > 1],100))+'\n')

clsizes_notone = [v for v in clusterlist.values() if v > 1]

fig = plt.figure()
minval = np.amin(clsizes_notone)
maxval = np.amax(clsizes_notone)
binstep = 100#(maxval-minval)/100
bins_ = np.arange(0,65001,binstep)#np.arange(0,maxval+0.000001,binstep)
hist, b = np.histogram(clsizes_notone,bins=bins_)
#hist = [float(i) for i in hist]
#for i in range(len(hist)):
#    if hist[i] > 0:
#        hist[i] = math.log10(hist[i])
x = bins_[:-1]+0.5*binstep
with open('hist_binsizes_max65000.txt','w') as f:
    for i in range(len(x)):
        f.write(str(x[i])+'\t'+str(hist[i])+'\n')
plt.plot(x,hist,color='b',markersize=1)
plt.savefig('hist_binsizes.png')
plt.close(fig)

'''
clsizes_notone_sorted = np.sort(clsizes_notone)
print(clsizes_notone_sorted[-10:])
print(clsizes_notone_sorted[-100])
print(clsizes_notone_sorted[-200])
print(clsizes_notone_sorted[-300])
print(clsizes_notone_sorted[-400])
print(clsizes_notone_sorted[-500])
print(clsizes_notone_sorted[-600])
print(clsizes_notone_sorted[-700])
print(clsizes_notone_sorted[-800])
print(clsizes_notone_sorted[-900])
print(clsizes_notone_sorted[-1000])
print(clsizes_notone_sorted[-1250])
print(clsizes_notone_sorted[-1500])
print(clsizes_notone_sorted[-1750])
'''
#clusterlist = None
'''
#----------------------
# Evaluate by scaffold
#----------------------

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom

read_scaffold_dict = json.load(open('/linuxhome/tmp/marleen/CAMIdataset1/read_scaffold_dict.json','r'))

tp = 0
fp = 0
no_scaffolds_per_cluster = {}
cluster_scaffolds = {}
for cl, reads in clusters.items():
    scaffolds = [read_scaffold_dict[rd] for rd in reads]
    un_scaffolds, count_scaffolds = np.unique(scaffolds,return_counts=True)
    no_scaffolds_per_cluster[cl] = len(un_scaffolds)
    cluster_scaffolds[cl] = []
    for i in range(len(un_scaffolds)):
        cluster_scaffolds[cl].append([un_scaffolds[i],count_scaffolds[i]])
    for c in count_scaffolds:
        tp += ncr(c,2)
    if len(un_scaffolds)>1:
        for i in range(len(count_scaffolds)):
            for j in range(i+1,len(count_scaffolds)):
                fp += count_scaffolds[i]*count_scaffolds[j]

no_clusters_per_scaffold = {}
scaffolds = np.unique([v for v in read_scaffold_dict.values()])
scaffold_clusters = {}
fn = 0
for scaf in scaffolds:
    scaffold_clusters[scaf] = []

for k,v in cluster_scaffolds.items():
    for sc in v:
        scaffold_clusters[sc[0]].append([k,sc[1]])

for scaf in scaffolds:
    no_clusters_per_scaffold[scaf] = len(scaffold_clusters[scaf])
    if no_clusters_per_scaffold > 1:
        for i in range(len(scaffold_clusters[scaf])):
            for j in range(i+1,len(scaffold_clusters[scaf])):
                fn += scaffold_clusters[scaf][i][1]*scaffold_clusters[scaf][j][1]

with open('cluster_evaluation_file_'+dataset+'.txt','w') as f:
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('Evaluation by scaffold\n')
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('Number of clusters per scaffold:\n')
    f.write('Min: '+str(np.amin([v for v in no_clusters_per_scaffold.values()]))+'\n')
    f.write('Mean: '+str(np.mean([v for v in no_clusters_per_scaffold.values()]))+'\n')
    f.write('Max: '+str(np.amax([v for v in no_clusters_per_scaffold.values()]))+'\n')
    f.write('\n')
    f.write('Number of scaffolds per cluster:\n')
    f.write('Min: '+str(np.amin([v for v in no_scaffolds_per_cluster.values()]))+'\n')
    f.write('Mean: '+str(np.mean([v for v in no_scaffolds_per_cluster.values()]))+'\n')
    f.write('Max: '+str(np.amax([v for v in no_scaffolds_per_cluster.values()]))+'\n')
    f.write('\n')
    f.write('Confusion matrix:\n')
    f.write('tp: '+str(tp)+'\n')
    f.write('fp: '+str(fp)+'\n')
    f.write('fn: '+str(fn)+'\n')
    f.write('\n')

read_scaffold_dict = None
no_scaffolds_per_cluster = None
cluster_scaffolds = None
no_clusters_per_scaffold = None
scaffolds = None
scaffold_clusters = None
tp = None
fp = None
fn = None


#----------------------
# Evaluate by strain
#----------------------

read_strain_dict = json.load(open('/linuxhome/tmp/marleen/CAMIdataset1/read_strain_dict.json','r'))

tp = 0
fp = 0
no_strains_per_cluster = {}
cluster_strains = {}
for cl, reads in clusters.items():
    strains = [read_strain_dict[rd] for rd in reads]
    un_strains, count_strains = np.unique(strains,return_counts=True)
    no_strains_per_cluster[cl] = len(un_strains)
    cluster_strains[cl] = []
    for i in range(len(un_strains)):
        cluster_strains[cl].append([un_strains[i],count_strains[i]])
    for c in count_strains:
        tp += ncr(c,2)
    if len(un_strains)>1:
        for i in range(len(count_strains)):
            for j in range(i+1,len(count_strains)):
                fp += count_strains[i]*count_strains[j]

no_clusters_per_strain = {}
strains = np.unique([v for v in read_strain_dict.values()])
strain_clusters = {}
fn = 0
for str in strains:
    strain_clusters[str] = []

for k,v in cluster_strains.items():
    for st in v:
        strain_clusters[st[0]].append([k,st[1]])

for str in strains:
    no_clusters_per_strain[str] = len(strain_clusters[str])
    if no_clusters_per_strain > 1:
        for i in range(len(strain_clusters[str])):
            for j in range(i+1,len(strain_clusters[str])):
                fn += strain_clusters[str][i][1]*strain_clusters[str][j][1]

with open('cluster_evaluation_file_'+dataset+'.txt','w') as f:
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('Evaluation by strain\n')
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('Number of clusters per strain:\n')
    f.write('Min: '+str(np.amin([v for v in no_clusters_per_strain.values()]))+'\n')
    f.write('Mean: '+str(np.mean([v for v in no_clusters_per_strain.values()]))+'\n')
    f.write('Max: '+str(np.amax([v for v in no_clusters_per_strain.values()]))+'\n')
    f.write('\n')
    f.write('Number of strains per cluster:\n')
    f.write('Min: '+str(np.amin([v for v in no_strains_per_cluster.values()]))+'\n')
    f.write('Mean: '+str(np.mean([v for v in no_strains_per_cluster.values()]))+'\n')
    f.write('Max: '+str(np.amax([v for v in no_strains_per_cluster.values()]))+'\n')
    f.write('\n')
    f.write('Confusion matrix:\n')
    f.write('tp: '+str(tp)+'\n')
    f.write('fp: '+str(fp)+'\n')
    f.write('fn: '+str(fn)+'\n')
    f.write('\n')

read_strain_dict = None
no_strains_per_cluster = None
cluster_strains = None
no_clusters_per_strain = None
strains = None
strain_clusters = None
tp = None
fp = None
fn = None
'''

#----------------------
# Evaluate by species
#----------------------

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom

read_species_dict = json.load(open(datafolder+'read_species_dict_S001.json','r'))

tp = 0
fp = 0
no_species_per_cluster = {}
cluster_species = {}
for cl, reads in clusters.items():
    species = [read_species_dict[rd] for rd in reads]
    un_species, count_species = np.unique(species,return_counts=True)
    no_species_per_cluster[cl] = len(un_species)
    cluster_species[cl] = []
    for i in range(len(un_species)):
        cluster_species[cl].append([un_species[i],count_species[i]])
    for c in count_species:
        tp += ncr(c,2)
    if len(un_species)>1:
        for i in range(len(count_species)):
            for j in range(i+1,len(count_species)):
                fp += count_species[i]*count_species[j]

with open('Clusters_with_multiple_species_toy1.txt','w') as f:
    for k,v in no_species_per_cluster.items():
        if v > 1:
            f.write('Cluster '+k+'\n')
            for sp in cluster_species[k]:
                f.write(str(sp)+'\n')

no_clusters_per_species = {}
species = np.unique([v for v in read_species_dict.values()])
species_clusters = {}
fn = 0
for spec in species:
    species_clusters[spec] = []

for k,v in cluster_species.items():
    for sp in v:
        species_clusters[sp[0]].append([k,sp[1]])

for spec in species:
    no_clusters_per_species[spec] = len(species_clusters[spec])
    if no_clusters_per_species[spec] > 1:
        for i in range(len(species_clusters[spec])):
            for j in range(i+1,len(species_clusters[spec])):
                fn += species_clusters[spec][i][1]*species_clusters[spec][j][1]

no_species_per_cluster_values = [v for v in no_species_per_cluster.values()]
no_clusters_per_species_values = [v for v in no_clusters_per_species.values()]

with open('cluster_evaluation_file_'+dataset+'.txt','a') as f:
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('Evaluation by species\n')
    f.write('-------------------------------------------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('Number of clusters per species:\n')
    f.write('Min: '+str(np.amin(no_clusters_per_species_values))+'\n')
    f.write('Mean: '+str(np.mean(no_clusters_per_species_values))+'\n')
    f.write('Max: '+str(np.amax(no_clusters_per_species_values))+'\n')
    f.write('0th percentile: '+str(np.percentile(no_clusters_per_species_values,0))+'\n')
    f.write('10th percentile: '+str(np.percentile(no_clusters_per_species_values,10))+'\n')
    f.write('20th percentile: '+str(np.percentile(no_clusters_per_species_values,20))+'\n')
    f.write('30th percentile: '+str(np.percentile(no_clusters_per_species_values,30))+'\n')
    f.write('40th percentile: '+str(np.percentile(no_clusters_per_species_values,40))+'\n')
    f.write('50th percentile: '+str(np.percentile(no_clusters_per_species_values,50))+'\n')
    f.write('60th percentile: '+str(np.percentile(no_clusters_per_species_values,60))+'\n')
    f.write('70th percentile: '+str(np.percentile(no_clusters_per_species_values,70))+'\n')
    f.write('80th percentile: '+str(np.percentile(no_clusters_per_species_values,80))+'\n')
    f.write('90th percentile: '+str(np.percentile(no_clusters_per_species_values,80))+'\n')
    f.write('100th percentile: '+str(np.percentile(no_clusters_per_species_values,100))+'\n')
    f.write('\n')
    f.write('Number of species per cluster:\n')
    f.write('Min: '+str(np.amin(no_species_per_cluster_values))+'\n')
    f.write('Mean: '+str(np.mean(no_species_per_cluster_values))+'\n')
    f.write('Max: '+str(np.amax(no_species_per_cluster_values))+'\n')
    f.write('0th percentile: '+str(np.percentile(no_species_per_cluster_values,0))+'\n')
    f.write('10th percentile: '+str(np.percentile(no_species_per_cluster_values,10))+'\n')
    f.write('20th percentile: '+str(np.percentile(no_species_per_cluster_values,20))+'\n')
    f.write('30th percentile: '+str(np.percentile(no_species_per_cluster_values,30))+'\n')
    f.write('40th percentile: '+str(np.percentile(no_species_per_cluster_values,40))+'\n')
    f.write('50th percentile: '+str(np.percentile(no_species_per_cluster_values,50))+'\n')
    f.write('60th percentile: '+str(np.percentile(no_species_per_cluster_values,60))+'\n')
    f.write('70th percentile: '+str(np.percentile(no_species_per_cluster_values,70))+'\n')
    f.write('80th percentile: '+str(np.percentile(no_species_per_cluster_values,80))+'\n')
    f.write('90th percentile: '+str(np.percentile(no_species_per_cluster_values,90))+'\n')
    f.write('100th percentile: '+str(np.percentile(no_species_per_cluster_values,100))+'\n')
    f.write('Fraction of clusters with more than 1 species: '+str(len([v for v in no_species_per_cluster_values if v > 1])/len(no_species_per_cluster_values))+'\n')
    f.write('Number of clusters with more than 1 species: '+str(len([v for v in no_species_per_cluster_values if v > 1]))+'\n')
    f.write('Number of clusters: '+str(len(no_species_per_cluster_values))+'\n')
    f.write('\n')
    f.write('Confusion matrix:\n')
    f.write('tp: '+str(tp)+'\n')
    f.write('fp: '+str(fp)+'\n')
    f.write('fn: '+str(fn)+'\n')
    f.write('\n')

with open('hist_clusters_per_species_'+dataset+'.txt','w') as f:
    minval = np.amin(no_clusters_per_species_values)
    maxval = np.amax(no_clusters_per_species_values)
    binstep = (maxval+1-minval)/100
    bins_ = np.logspace(np.log10(1),np.log10(maxval),50)
#    bins_ = np.arange(0,maxval+1,binstep)
    no_clusters_per_species_values = [float(i) for i in no_clusters_per_species_values]
    hist, b = np.histogram(no_clusters_per_species_values,bins=bins_)
    hist = [float(i) for i in hist]
    x = bins_[:-1]+0.5*binstep
    for i in range(len(x)):
        f.write(str(x[i])+'\t'+str(hist[i])+'\n')

with open('list_clusters_per_species_'+dataset+'.txt','w') as f:
    for k,v in no_clusters_per_species.items():
        f.write('Species '+str(k)+'\t'+str(v)+' clusters\n')

with open('hist_species_per_cluster_'+dataset+'.txt','w') as f:
    minval = np.amin(no_species_per_cluster_values)
    maxval = np.amax(no_species_per_cluster_values)
    binstep = 1 #(maxval+1-minval)/100
    bins_ = np.arange(0.5,maxval+1,binstep)
    hist, b = np.histogram(no_species_per_cluster_values,bins=bins_)
    hist = [float(i) for i in hist]
    x = bins_[:-1]+0.5*binstep
    for i in range(len(x)):
        f.write(str(x[i])+'\t'+str(hist[i])+'\n')


with open('cluster_evaluation_file_'+dataset+'.txt','a') as f:
    with open('data_for_excel_'+dataset+'.txt','w') as fm:
        with open('data_for_coverage_hist_'+dataset+'.txt','w') as fh:
            with open('data_for_clustsperspecies_bar_'+dataset+'.txt','w') as fb:
                f.write('----------------------------------------------------------------------\n')
                f.write('Min, mean and max clustersizes per species\n')
                f.write('----------------------------------------------------------------------\n')
                fm.write('\t\tNo. reads of this species in cluster\n')
                fm.write('Species\tNumber of read pairs\tMin\tMean\tMax\n')
                for k,v in species_clusters.items():
                    sizelist = []
                    for cl in v:
                        sizelist.append(cl[1])
                    if len(sizelist) > 0:
                        f.write('Species '+k+' has '+str(len(sizelist))+' clusters with no. read pairs: min '+str(np.amin(sizelist))+', mean '+str(np.mean(sizelist))+', max '+str(np.amax(sizelist))+', total '+str(sum(sizelist))+'. \n')
                        fm.write(k+'\t'+str(sum(sizelist))+'\t'+str(np.amin(sizelist))+'\t'+str(np.mean(sizelist))+'\t'+str(np.amax(sizelist))+'\n')
                        fh.write(k+'\t'+str(sum(sizelist))+'\n')
                        fb.write(k+'\t'+str(len(sizelist))+'\n')
                    else:
                        f.write('Species '+k+': no clusters. \n')
                        fm.write(k+'\t0\t0\t0\t0\t0\n')
                        fh.write(k+'\t0\n')
                        fb.write(k+'\t0\n')


reads_per_species = json.load(open(datafolder+'num_reads_per_species_S001.json','r'))
reads_clustered = {}
frac_reads_clustered = {}
for k,v in species_clusters.items():
    sizelist = []
    for cl in v:
        sizelist.append(cl[1])
    reads_clustered[k] = sum(sizelist)
    frac_reads_clustered[k] = float(reads_clustered[k])/float(reads_per_species[k])

json.dump(frac_reads_clustered,open('frac_reads_clustered_'+dataset+'.json','w'))

'''
covs = []
fracs = []
for k in frac_reads_clustered.keys():
    covs.append(sp_cov[k])
    fracs.append(frac_reads_clustered[k])

fig = plt.figure()
plt.scatter(covs,fracs,s=1)
xmin = math.floor(np.amin(covs))
xmax = math.ceil(np.amax(covs))
ymin = math.floor(np.amin(fracs))
ymax = math.ceil(np.amax(fracs))
plt.xticks(np.arange(xmin,xmax+1,10))
plt.yticks(np.arange(ymin,ymax+0.01,0.1))
plt.xlabel('Coverage')
plt.ylabel('Fraction reads clustered')
plt.savefig('covplot_'+dataset+'.png')
plt.close(fig)
'''
