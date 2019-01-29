import json
import sys

fqfile = sys.argv[1]
skip_digits = int(sys.argv[2])+1
subset = sys.argv[3]

read_seq_phred = {}
read_len = {}

with open(fqfile,'r') as f:
    i = 0
    for line in f:
        if i % 4 == 0:
            rd = line.rstrip()[skip_digits:]
        elif i % 4 ==1:
            seq = line.rstrip()
        elif i % 4 == 3:
            phred = line.rstrip()
            read_seq_phred[rd] = [seq,phred]
            read_len[rd] = len(seq)
            rd = None
            seq = None
            phred = None
        i += 1

json.dump(read_seq_phred,open('read_seq_phred_'+subset+'.json','w'))
json.dump(read_len,open('read_len_dict_'+subset+'.json','w'))
