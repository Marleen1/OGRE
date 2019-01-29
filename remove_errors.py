import json
import sys
import os
import multiprocessing as mp
import math

folder = sys.argv[1]
datafolder = sys.argv[2]
subset = sys.argv[3]

read_species = json.load(open(datafolder+'read_species_dict_'+subset+'.json','r'))

with open(folder+'ovlp.txt','r') as fr:
    with open(folder+'ovlp_corrected.txt','w') as fw:
        with open(folder+'ovlp_removedlines.txt','w') as fe:
            for line in fr:
                ln = line.rstrip().split('\t')
                if ln[0][:-2] in read_species.keys() and ln[1][:-2] in read_species.keys():
                    fw.write(line)
                else:
                    fe.write(line)
