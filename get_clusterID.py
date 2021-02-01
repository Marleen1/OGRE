import json
import sys

cluster_file,outdir=sys.argv[1:]

d={}

with open(cluster_file) as fr:
    data=json.load(fr)
    for key,val in data.items():
        val2=[v.split("|")[-1] for v in val]
        d[key]=val2

with open(outdir+"/cluster2reads.json",'w') as fw:
    json.dump(d,fw)

with open(outdir+"/cluster.id",'w') as fw:
    fw.write('\n'.join(d.keys()))

