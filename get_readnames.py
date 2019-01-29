import sys

folder = sys.argv[1]
file = sys.argv[2]

i = 0
with open(folder+file,'r') as fr:
    with open('readnames.txt','w') as fw:
        for line in fr:
            if i % 8 == 0 :
                fw.write(line[1:-3]+'\n')
            i += 1
