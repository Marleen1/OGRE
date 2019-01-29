import math

numsplits = 60
folder = ".../CAMIdata/CAMI3/"
fqfile = "RH_S001__insert_270_qtrimmed"
numlines = 50000000
lines_per_file = 2500000
subset = "S001"
train_data = 'traindata_CAMI12toy23_80000.txt'
val_data = 'CAMI3'
run_ID = 'CAMI3_sr_k21_w11_s60_m60_n2_r0_A4_B2'

num_splits = math.ceil(numlines/lines_per_file)
numsplitlines = math.ceil(numlines/(8*numsplits))*8
skip_digits = '0'


files = []
for g1 in range(1,numsplits+1):
    for g2 in range(g1,numsplits+1):
        files.append([str(g1).zfill(2),str(g2).zfill(2)])

singles = []
for g1 in range(1,numsplits+1):
    singles.append(str(g1).zfill(2))

splits = []
for g1 in range(1,num_splits+1):
    splits.append(str(g1).zfill(4))


rule all:
    input:
        'covplot_'+run_ID+'.png'

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
        shell("python make_read_seq_phred.py {input} "+skip_digits+" g{wildcards.f}")


rule minimap:
    input:
        folder+fqfile+"_g{f0}.fq",
        folder+fqfile+"_g{f1}.fq",
        "read_seq_phred_g{f0}.json",
        "read_seq_phred_g{f1}.json"
    output:
        "ovlp_g{f0}_g{f1}_summ.txt"
    threads:
        4
    run:
        shell("minimap2 -t {threads} --sr -X -a -k 21 -w 11 -s 60 -m 60 -n 2 -r 0 -A 4 -B 2 --end-bonus=100 {input[0]} {input[1]} > ovlp_g{wildcards.f0}_g{wildcards.f1}.sam")
        shell("python skip_sam_header.py ovlp_g{wildcards.f0}_g{wildcards.f1}.sam ovlp_g{wildcards.f0}_g{wildcards.f1}_nohead.sam")
        shell("/bin/rm ovlp_g{wildcards.f0}_g{wildcards.f1}.sam")
        shell("python keep_relevant_data_sam.py ovlp_g{wildcards.f0}_g{wildcards.f1}_nohead.sam {threads} "+skip_digits+" {output} "+folder+" ./ g{wildcards.f0} g{wildcards.f1}"),
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
        4
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

rule split_same_diff_spec:
    input:
        "ovlp.txt"
    output:
        "ovlp_samespec.txt",
        "ovlp_diffspec.txt"
    threads:
        1
    run:
        shell("python remove_errors.py ./ "+folder+" "+subset),
        shell("/bin/rm ovlp.txt"),
        shell("mv ovlp_corrected.txt ovlp.txt"),
        shell("python same_or_diff_species.py ./ "+folder+" {threads} "+subset)

rule fitModel:
    input:
        train_data
    output:
        'ovlp_LogRegr_fittedModel_'+val_data+'.sav'
    threads:
        18
    run:
        shell('python logregr.py {input} {output} {threads}')

rule split2:
    input:
        'ovlp_LogRegr_fittedModel_'+val_data+'.sav',
        'ovlp.txt'
    output:
        expand('ovlp_split{s}.txt', s=splits)
    threads:
        18
    run:
        shell('split -l '+str(lines_per_file)+' -a 4 --numeric-suffixes=1 --additional-suffix=.txt ovlp.txt ovlp_split')


rule predict:
    input:
        'ovlp_split{s}.txt',
        'ovlp_LogRegr_fittedModel_'+val_data+'.sav'
    output:
        'ovlp_predict_keep_split{s}.txt',
        'ovlp_predict_removed_split{s}.txt'
    threads:
        1
    run:
        shell('python predict_samediffspec.py {input[1]} {wildcards.s}'),
        shell('/bin/rm ovlp_split{wildcards.s}.txt')

rule merge_keep:
    input:
        expand('ovlp_predict_keep_split{s}.txt', s=splits)
    output:
        'ovlp_predictLogRegr_'+val_data+'.txt'
    threads:
        18
    run:
        shell('cat ovlp_predict_keep_split* > ovlp_predict_temp.txt')
        shell('/bin/rm ovlp_predict_keep_split*')
        shell('sort -nr -k3,3 -T ./temp ovlp_predict_temp.txt > {output}')
        shell('/bin/rm ovlp_predict_temp.txt')

rule readnames:
    input:
        folder+fqfile+".fq"
    output:
        'readnames.txt'
    threads:
        1
    run:
        shell('python get_readnames.py '+folder+' '+fqfile+'.fq')

rule clustering:
    input:
        'ovlp_predictLogRegr_'+val_data+'.txt',
        'readnames.txt'
    output:
        'cluster_evaluation_file_'+run_ID+'_max{limit}_final.txt'
    threads:
        36
    run:
        shell('python bin_pointer_limited_filechunks_shortpath.py {input} readnames.txt {wildcards.limit} '+run_ID+' {threads}'),
        shell('/bin/rm Chunkfile*'),
        shell('python getclusters.py '+run_ID+'_max{wildcards.limit}_final {threads}'),
        shell('python evaluate_clusters.py '+run_ID+'_max{wildcards.limit}_final '+run_ID+'_max{wildcards.limit}_final_clustersizes.json '+folder)

rule combined_output:
    input:
        expand('cluster_evaluation_file_'+run_ID+'_max{limit}_final.txt', limit=['inf'])
    output:
        'covplot_'+run_ID+'.png'
    threads:
        1
    run:
        shell('python covplot.py '+run_ID+' '+folder)

