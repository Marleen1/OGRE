import json
import sys
import os
import multiprocessing as mp
import math

folder = sys.argv[1]
datafolder = sys.argv[2]
threads = int(sys.argv[3])
subset = sys.argv[4]
infile = sys.argv[5]
#subset = ''

def splitPairTypes(fileID):
    with open('tempfile_split'+str(fileID).zfill(2)+'.paf','r') as fr:
        with open('tempfile_samespec'+str(fileID).zfill(2)+'.paf','w') as fsame:
            with open('tempfile_diffspec'+str(fileID).zfill(2)+'.paf','w') as fdiff:
                for line in fr:
                    ln = line.rstrip().split('\t')
                    if any(x in read_species[ln[0][:-2]] for x in read_species[ln[1][:-2]]):
                        fsame.write(line)
                    else:
                        fdiff.write(line)
    return


if __name__ == '__main__':
#    read_species = json.load(open(datafolder+'read_species_dict.json','r'))
    read_species = json.load(open(datafolder+'read_species_dict_'+subset+'.json','r'))
    with open(folder+infile,'r') as f:
        num_lines = sum(1 for line in f)
    num_lines_split = math.ceil(num_lines/threads)
    os.system('split -l '+str(num_lines_split)+' -a 2 --numeric-suffixes=1 --additional-suffix=.paf '+folder+infile+' tempfile_split')
    po = mp.Pool(threads)
    for th in range(1,threads+1):
        po.apply_async(splitPairTypes, args=(th,))
    po.close()
    po.join()
    os.system('cat tempfile_samespec* > '+folder+'ovlp_predictLogRegr_samespec.txt')
    os.system('cat tempfile_diffspec* > '+folder+'ovlp_predictLogRegr_diffspec.txt')
    os.system('/bin/rm tempfile_split*')
    os.system('/bin/rm tempfile_samespec*')
    os.system('/bin/rm tempfile_diffspec*')
