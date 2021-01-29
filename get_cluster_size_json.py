#!/usr/bin/env python

import json
import pysam
import pandas as pd
import os
import random
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import sys


cluster_json=sys.argv[1]
outfile=sys.argv[2]

cluster2read=json.load(open(cluster_json,'r'))

cluster2size= {cluster:len(reads) for cluster,reads in cluster2read.items()}

json.dump(cluster2size,open(outfile,'w'))

