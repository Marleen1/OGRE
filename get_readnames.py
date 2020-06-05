import sys
import numpy as np

folder = sys.argv[1]
fqfile = sys.argv[2]
#runID = sys.argv[3]
#skip_digits = int(sys.argv[3])
#subset = sys.argv[4]

with open(folder+fqfile,'r') as fr:
    with open('readnames.txt','w') as fw:
        i = 0
        for line in fr:
            if i % 8 == 0:
                fw.write(line.rstrip().replace('@','')[:-2]+'\n')
            i += 1
