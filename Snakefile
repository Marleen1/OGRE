import math

folder = "/linuxhome/tmp/marleen/CAMI/"
fqfile = "RL_S001__insert_270_qtrimmed_spec340_223_220_126.fq"
limits = ['inf'] #,'3300','17000','33000']

fqfile = fqfile.replace('.fq','')


train_data = 'traindata_CAMI23toy23_80000.txt'
with open(folder+fqfile+".fq","r") as f:
    numlines = sum(1 for line in f)
numsplitlines = 8000000 
numsplits = math.ceil(numlines/numsplitlines)
lines_per_file = 2500000

files = []
for g1 in range(1,numsplits+1):
    for g2 in range(g1,numsplits+1):
        files.append([str(g1).zfill(2),str(g2).zfill(2)])

singles = []
for g1 in range(1,numsplits+1):
    singles.append(str(g1).zfill(2))

rule all:
    input:
        expand(fqfile+'_max{limit}_final_clusters_grouped.json',limit=limits)

rule split:
    input:
        folder+fqfile+".fq"
    output:
        expand(folder+fqfile+"_g{f}.fq", f=singles)
    run:
        shell("split -l "+str(numsplitlines)+" --numeric-suffixes=1 --additional-suffix=.fq {input} "+folder+fqfile+"_g")


rule read_seq_phred:
    input:
        folder+fqfile+"_g{f}.fq"
    output:
        "read_seq_phred_g{f}.json"
    threads:
        1
    run:
        shell("python make_read_seq_phred.py {input} g{wildcards.f}")


rule minimap:
    input:
        folder+fqfile+"_g{f0}.fq",
        folder+fqfile+"_g{f1}.fq",
        "read_seq_phred_g{f0}.json",
        "read_seq_phred_g{f1}.json"
    output:
        "ovlp_g{f0}_g{f1}_summ.txt"
    threads:
        12
    run:
        shell("minimap2 -t {threads} --sr -X -a -k 21 -w 11 -s 60 -m 60 -n 2 -r 0 -A 2 -B 2 --end-bonus=100 {input[0]} {input[1]} > ovlp_g{wildcards.f0}_g{wildcards.f1}.sam")
        shell("python skip_sam_header.py ovlp_g{wildcards.f0}_g{wildcards.f1}.sam ovlp_g{wildcards.f0}_g{wildcards.f1}_nohead.sam")
        shell("/bin/rm ovlp_g{wildcards.f0}_g{wildcards.f1}.sam")
        shell("python keep_relevant_data_sam.py ovlp_g{wildcards.f0}_g{wildcards.f1}_nohead.sam {threads} {output} "+folder+" ./ g{wildcards.f0} g{wildcards.f1}"),
        shell("/bin/rm ovlp_g{wildcards.f0}_g{wildcards.f1}_nohead.sam")


rule combine:
    input:
        "ovlp_g{f0}_g{f1}_summ.txt",
        "ovlp_g{f1}_g{f0}_summ.txt"
    output:
        "ovlp_g{f0}_g{f1}_c.txt"
    threads:
        1
    run:
        if wildcards.f0 == wildcards.f1:
            shell("cat {input[0]} > {output}"),
            shell("/bin/rm {input[0]}")
        else:
            shell("cat {input} > {output}"),
            shell("/bin/rm {input}")

rule nodoubles:
    input:
        "ovlp_g{f0}_g{f1}_c.txt"
    output:
        "ovlp_g{f0}_g{f1}_combined.txt"
    threads:
        8
    run:
        shell("awk '!seen[$0]++' {input} > {output}")
        shell("/bin/rm {input}")

rule concatall:
    input:
        expand("ovlp_g{f[0]}_g{f[1]}_combined.txt", f=files)
    output:
        "ovlp.txt"
    threads:
        8
    run:
        shell("cat {input} > {output}")
        shell("/bin/rm {input}")


rule fitModel:
    input:
        train_data
    output:
        'ovlp_LogRegr_fittedModel.sav'
    threads:
        12
    run:
        shell('python logregr.py {input} {output} {threads}')

checkpoint split_2:
    input:
        'ovlp.txt'
    output:
        splitd = directory('ovlp_split')
    run:
        shell('mkdir {output.splitd}'),
        shell('split -l '+str(lines_per_file)+' -a 4 --numeric-suffixes=1 --additional-suffix=.txt ovlp.txt {output.splitd}/')

rule predict:
    input:
        'ovlp_split/{s}.txt',
        'ovlp_LogRegr_fittedModel.sav'
    output:
        'ovlp_predict_split{s}.txt',
    threads:
        2
    run:
        shell('python predict_samediffspec.py {input[1]} {wildcards.s}'),
        #shell('/bin/rm {input[0]}')

def aggregate_split(wildcards):
    chkp_done = checkpoints.split_2.get().output.splitd
    chkp_output = sorted(glob_wildcards(os.path.join(chkp_done, "{s}.txt")).s)
    return expand('ovlp_predict_split{s}.txt', s= chkp_output)

rule merge_keep:
    input:
        aggregate_split
    output:
        'ovlp_predictLogRegr.txt'
    threads:
        36
    run:
        shell('mkdir temp')
        shell('cat ovlp_predict_split* > ovlp_predict_temp.txt')
        shell('/bin/rm ovlp_predict_split*')
        shell('sort -nr -k3,3 -T ./temp ovlp_predict_temp.txt > {output}')
        shell('/bin/rm -r temp')
        shell('/bin/rm ovlp_predict_temp.txt')


rule get_readnames:
    output:
        'readnames.txt'
    threads:
        1
    run:
        shell('python get_readnames.py '+folder+' '+fqfile+'.fq')

rule clustering:
    input:
        'ovlp_predictLogRegr.txt',
        'readnames.txt'
    output:
        fqfile+'_max{limit}_final_clusters_grouped.json'
    threads:
        72
    run:
        shell('python bin_pointer_limited_filechunks_shortpath.py {input[0]} {wildcards.limit} '+fqfile+' {threads}'),
        shell('/bin/rm Chunkfile*'),
        shell('python getclusters.py '+fqfile+'_max{wildcards.limit}_final {threads}'),
