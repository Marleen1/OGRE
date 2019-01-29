import numpy as np
import csv
import json
import time
import sys
import os
import multiprocessing
import math
from multiprocessing import Value
import collections
import os
from pathlib import Path

#print 'Run this script in Python 2'

threads = int(sys.argv[5])
#os.remove('Chunkfile')
starttime = time.time()

ovlpfile = sys.argv[1] #'/linuxhome/tmp/marleen/Minimap/ovlp_predictLogRegr.txt' #'ovlp_sr_s40_w8_k12_m40_n2_g01_g01_combined_sorted_short55.paf' #'testovlp'
readnamefile = sys.argv[2] #'/linuxhome/tmp/marleen/CAMI/CAMIdataset3/readnames_S001.txt' #'readnames_g01_g01_sorted' #'testreadnames'
maxsize = sys.argv[3] #33000 #float('inf')
rID = sys.argv[4]
run_ID = rID+'_max'+maxsize
if maxsize == 'inf':
    maxsize = float('inf')
else:
    maxsize = int(maxsize)
#num_lines_processed_together = 100000
#chunksize = num_files_processed_together/threads
chunksize = 2600000 #This is approximately the number of bytes used for 100000 lines

def chunkify(fname):
    fileEnd = os.path.getsize(fname)
    with open(fname,'rb') as f: #must be binary rb
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(chunksize,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def findhead_nochanges(r):
    stop = 0
    rcurrent = r
    while stop == 0:
        if isinstance(clusters[rcurrent], int):
            clust = clusters[rcurrent]
            stop = 1
        else:
            rcurrent = clusters[rcurrent]
    return (rcurrent,clust)

def findhead_lim(r):
    stop = 0
    rcurrent = r
    pathlen = 1
    adapt = []
    while stop == 0:
        if isinstance(clusters[rcurrent], int):
            clust = clusters[rcurrent]
#            if pathlen >= 3:
#                for node in adapt[:-1]:
#                    clusters[node] = rcurrent
            stop = 1
        else:
            adapt.append(rcurrent)
            rcurrent = clusters[rcurrent]
            pathlen += 1
    return (rcurrent,clust,pathlen)

def clusteralgorithm(fileID):#(fileID,clusters,clusterlist):
    linecount = 0
    with open(fileID,'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            matches = row[1]
            cllim_r1_head, clust1, pathlen1 = findhead_lim(row[0][:-2])
            cllim_r2_head, clust2, pathlen2 = findhead_lim(row[1][:-2])
            if clust1 != clust2:
                clustsize = clusterlist[clust1]+clusterlist[clust2]
                if clustsize <= maxsize:
                    if pathlen2 < pathlen1:
                        clusters[cllim_r2_head] = cllim_r1_head
                        clusterlist[clust1] = clustsize
                        clusterlist.pop(clust2)
                    else:
                        clusters[cllim_r1_head] = cllim_r2_head
                        clusterlist[clust2] = clustsize
                        clusterlist.pop(clust1)
            linecount += 1
    return linecount, matches

def getchunkfile(a):#(chunknum_,chunkstart_,chunksize_):
    chunknum_ = a[0]
    chunkstart_ = a[1]
    chunksize_ = a[2]
    with open('Chunkfile_'+run_ID+'_'+str(chunknum_),'w') as fw:
        with open(ovlpfile,'r') as fo: #Must be binary?
            fo.seek(chunkstart_)
            lines = fo.read(chunksize_).splitlines()
            for line in lines:
                linesplit = line.rstrip().split('\t')
                cllim_r1_head, clust1 = findhead_nochanges(linesplit[0][:-2])
                cllim_r2_head, clust2 = findhead_nochanges(linesplit[1][:-2])
                if clust1 != clust2 and clusterlist[clust1]+clusterlist[clust2] < maxsize:
                    fw.write(line+'\n')
    return

nmatches = 150
if __name__ == '__main__':
    clusters = {}
    clusterlist = {}
    i = 1
    with open(readnamefile,'r') as f:
        for line in f:
            clusters[line.rstrip()] = i
            clusterlist[i] = 1
            i += 1
    # Make chunks and sessions
    starttime = time.time()
    sessiondict = {}
    session = 1
    chunknum = 1
    sessionlist = []
    for chunkStart,chunkSize in chunkify(ovlpfile):
        sessiondict[session]=[session,chunkStart,chunkSize]
        session += 1
#    print('Chunks made' )
    numsessions = len(sessiondict.keys())
    starttime = time.time()
    sess = 1
    while sess <= numsessions:
        po = multiprocessing.Pool(threads)
        for th in range(threads):
            nfiles = 0
            if sess <= numsessions:
                session = sessiondict[sess]
                sess += 1
                nfiles += 1
                po.apply_async(getchunkfile, args=([th+1,session[1],session[2]],))
        po.close()
        po.join()
        with open('Chunkfile_'+run_ID,'w') as fcat:
            for ch in range(1,threads+1):
                my_file = Path('Chunkfile_'+run_ID+'_'+str(ch))
                if my_file.exists():
                    with open('Chunkfile_'+run_ID+'_'+str(ch),'r') as infile:
                        fcat.write(infile.read())
        if os.stat('Chunkfile_'+run_ID).st_size != 0:
            linecount, matches = clusteralgorithm('Chunkfile_'+run_ID)
            if sess % threads == 100:
                with open(ovlpfile,'r') as fo:
                    fo.seek(sessiondict[sess][1])
                    for line in fo:
                        ln = line.split('\t')
                        threshold = ln[2].rstrip()
                        break
                with open(run_ID+'_round'+str(sess)+'_threshold'+str(threshold)+'_clusters.json','w') as fp:
                    json.dump(clusters,fp)
                with open(run_ID+'_round'+str(sess)+'_threshold'+str(threshold)+'_clustersizes.json','w') as fp:
                    json.dump(clusterlist,fp)
#            if int(matches) < int(nmatches):
#                nmatches = matches
#                with open('minimap_readclusterdict_until'+str(matches)+'_unlimited.json','w') as fp:
#                    json.dump(clusters, fp)
#                with open('minimap_readclustersizes_until'+str(matches)+'_unlimited.json','w') as fp:
#                    json.dump(clusterlist, fp)
        my_file = Path('Chunkfile_'+run_ID)
        if my_file.exists():
            os.remove('Chunkfile_'+run_ID)
        for ch in range(1,nfiles+1):
            my_file = Path('Chunkfile_'+run_ID+'_'+str(ch))
            if my_file.exists():
                os.remove('Chunkfile_'+run_ID+'_'+str(ch))
        nextline = (sess-1)*100000
#        print('About %f lines processed in %f seconds.' % (nextline,time.time()-starttime))
#        print('Chunkfile contained %f lines, so %f percent.' % (linecount,100*linecount/(100000*threads)))
#        print('%f out of %f sessions done.' % (sess-1,numsessions))
#        if sess <= numsessions:
#            with open(ovlpfile,'r') as fopen:
#                fopen.seek(sessiondict[sess][1])
#                print(fopen.readline().rstrip())

#    print('Number of clusters with size limit: '+str(len(clusterlist.keys())))
#    print('Min cluster size of clusters with size limit: '+str(np.amin([val for val in clusterlist.values()])))
#    print('Mean cluster size of clusters with size limit: '+str(np.mean([val for val in clusterlist.values()])))
#    print('Max cluster size of clusters with size limit: '+str(np.amax([val for val in clusterlist.values()])))
    with open(run_ID+'_final_clusters.json','w') as fp:
        json.dump(clusters,fp)
    with open(run_ID+'_final_clustersizes.json','w') as fp:
        json.dump(clusterlist,fp)

#         Run clusteralgorithm


'''
f = open('ovlp_sr_s40_w8_k12_m40_n2_all_sorted_viasteps_short40.paf','rb')
ct = time.time()
f.readline()
start = f.tell()
print(f.readline())
f.close()



    starttime = time.time()
    po = multiprocessing.Pool(threads)
    for matches in range(150,139,-1):
        print(matches)
        for g1 in range(1,41):
            gname = 'g'+str(g1).zfill(2)+'_g'+str(g1).zfill(2)
            po.apply_async(clusteralgorithm(gname,matches,clusters,clustersizes,int(maxsize/40)))
    po.close()
    po.join()
    print(time.time()-starttime)
    totclusts = 0
    for key in clustersizes.keys():
        totclusts += len(clustersizes[key])
    print(totclusts)
    print(time.time()-starttime)
    print(clusters['g01_g01']['R1160995'])
    print(clusters['g01_g01']['R694278'])
    print(clusters['g01_g01']['R1065371'])
    print(clusters['g01_g01']['R393621'])
'''
'''
with open('minimap_readclusterdict_until'+str(nmatches)+'_unlimited.json','w') as fp:
    json.dump(clusters_limited, fp)
with open('minimap_clusterreadsdict_until'+str(nmatches)+'_unlimited.json','w') as fp:
    json.dump(clusterlist_limited, fp)
'''
'''
contignames = list(set([key for key in clusterlist_limited.keys()]))

patches = {}
for contig in contignames:
    patches[contig] = [key for key, val in clusters_limited.items() if clusters_limited[key][0]==contig]

with open('minimap_clusternames_limited.txt','w') as f:
    for name in contignames:
        f.write('patch_'+str(name)+'\n')

with open('minimap_readclusterdict_limited.json', 'w') as fp:
    json.dump(clusters_limited, fp)

with open('minimap_patches_limited.json','w') as fp:
    json.dump(patches, fp)

contignames = list(set([key for key in clusters.keys()]))

patches = {}
for contig in contignames:
    patches[contig] = [key for key, val in clusters.items() if clusters[key][0]==contig]

with open('minimap_clusternames.txt','w') as f:
    for name in contignames:
        f.write('patch_'+str(name)+'\n')

with open('minimap_readclusterdict.json', 'w') as fp:
    json.dump(clusters, fp)

with open('minimap_patches.json','w') as fp:
    json.dump(patches, fp)
'''
'''
print('Number of clusters without size limit: '+str(len(clusterlist.keys())))
print('Min cluster size of clusters without size limit: '+str(np.amin(clustsizes)))
print('Mean cluster size of clusters without size limit: '+str(np.mean(clustsizes)))
print('Max cluster size of clusters without size limit: '+str(np.amax(clustsizes)))

print('Number of clusters with size limit: '+str(len(clusterlist_limited.keys())))
print('Min cluster size of clusters with size limit: '+str(np.amin(clustsizes_limited)))
print('Mean cluster size of clusters with size limit: '+str(np.mean(clustsizes_limited)))
print('Max cluster size of clusters with size limit: '+str(np.amax(clustsizes_limited)))
'''
