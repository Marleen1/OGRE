import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import json

dataset = sys.argv[1]
datafolder = sys.argv[2]

limits = ['inf']
labels = ['unlimited']
colors = ['r','b','g','k']

#covs = [10,20,32,55,63]
covs = json.load(open(datafolder+'species_cov.json','r'))
#fracses = [[0.2,0.3,0.4,0.5,0.6],[0.25,0.33,0.27,0.45,0.5],[0.22,0.32,0.42,0.53,0.64],[0.33,0.28,0.47,0.52,0.58]]
fig = plt.figure()
for i in range(len(limits)):
#    fracs = fracses[i]
    fracs = json.load(open('frac_reads_clustered_'+dataset+'_max'+limits[i]+'_final.json','r'))
    covlist = []
    fraclist = []
    for k,v in covs.items():
        covlist.append(v)
        fraclist.append(fracs[k])
    plt.scatter(covlist,fraclist,s=7,c=colors[i],label=labels[i])

xmin = math.floor(np.amin([v for v in covs.values()]))
xmax = math.ceil(np.amax([v for v in covs.values()]))
ymin = 0 #math.floor(np.amin(fracs.keys()))
ymax = 1 #math.ceil(np.amax(fracs.keys()))
plt.xticks(np.arange(xmin,xmax+1,(xmax-xmin)/5))
plt.yticks(np.arange(ymin,ymax+0.01,0.1))
plt.legend()
plt.xlabel('Coverage')
plt.ylabel('Fraction reads clustered')
plt.savefig('covplot_'+dataset+'.png')
plt.close(fig)

