#!/usr/bin/env python
import json
import re
import random
import sys,os


fastq = sys.argv[1]
cluster_json = sys.argv[2]
sub_file = sys.argv[3]
if_merge = sys.argv[4]
prefix=sys.argv[5]
outdir0=sys.argv[6]

tar_clusters=[]
with open(sub_file,'r') as fr:
    for line in fr:
        tar_clusters.append(line.strip())

cluster2read =json.load(open(cluster_json,'r'))
cluster2read = {c:reads for c,reads in cluster2read.items() if c in tar_clusters}


read2cluster = {read:key for key,val in cluster2read.items() for read in val}

cluster_len=len(tar_clusters)

def get_fq4cluster(tar_clusters,file,outdir):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["fq1w%s"%cluster] = open(outdir4c + cluster + '.1.fq', 'w')
        locals()["fq2w%s"%cluster] = open(outdir4c + cluster + '.2.fq', 'w')
    #
    pre_record = []
    pre_read = ''
    pre_header = ''
    if_paired = 0
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            #             if i>100000:break
            #if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 4 == 0:
                if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
                    if re.search(r'/1$', pre_header):
                        locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
                        if_paired = 0
                    else:
                        locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
                        if_paired = 1
                pre_record = []
                pre_record.append(line)
                pre_read = re.split('[@|/]', line)[-2]
                #pre_read = re.split('[@/]', line)[-2]
                pre_header = line
            else:
                pre_record.append(line)
        if not if_paired: locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
    fq.close()
    for cluster in tar_clusters:
        locals()["fq1w%s"%cluster].close()
        locals()["fq2w%s"%cluster].close()

def get_fq4cluster2(tar_clusters,file,outdir):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["fq1w%s"%cluster] = open(outdir4c + cluster + '.fq', 'w')
    #
    pre_record = []
    pre_read = ''
    pre_header = ''
    if_paired = 0
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            #             if i>100000:break
            #if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 4 == 0:
                if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
                    if re.search(r'/1$', pre_header):
                        locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
                        if_paired = 0
                    else:
                        locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
                        if_paired = 1
                pre_record = []
                pre_record.append(line)
                #pre_read = re.split('[@/]', line)[-2]
                pre_read = re.split('[@|/]', line)[-2]
                pre_header = line
            else:
                pre_record.append(line)
        if (not if_paired) and (pre_read in read2cluster):
            locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
    fq.close()
    for cluster in tar_clusters:
        locals()["fq1w%s"%cluster].close()


#print("begin...")
#print("#"*50)

#tar_clusters=sys.argv[1]
k=420 if if_merge != 'true' else 840
split_n = int(cluster_len/k)+1
outdir=outdir0+"/fastq_"+prefix
for i in range(split_n):
    #print("the "+str(i+1)+"/"+str(split_n)+" part start...")
    start=k*i
    end=k*(i+1) if (k*(i+1))<cluster_len else  cluster_len
    cc=tar_clusters[start:end]
    if if_merge != 'true':
        get_fq4cluster(cc,fastq,outdir)
    else:
        get_fq4cluster2(cc,fastq,outdir)
    #print("the "+str(i+1)+"/"+str(split_n)+" part finished...\n")
    #print("#"*50)

print("done...")

