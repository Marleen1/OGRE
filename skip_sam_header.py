import numpy as np
import sys

samin = sys.argv[1]
samout = sys.argv[2]

with open(samin,'r') as fr:
    with open(samout,'w') as fw:
        for line in fr:
            if line[0] != '@':
                ln = line.split('\t')
                if ln[2] != '*':
                    fw.write(line)
                break
        for line in fr:
            ln = line.split('\t')
            if ln[2] != '*':
                fw.write(line)
