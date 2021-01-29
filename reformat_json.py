import json
d={}
with open("./cluster2reads.json2") as fr:
    data=json.load(fr)
    for key,val in data.items():
        val2=[v.split("|")[-1] for v in val]
        d[key]=val2
with open("./cluster2reads.json",'w') as fw:
    json.dump(d,fw)

