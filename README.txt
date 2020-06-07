OGRE is a read clustering tool that clusters short reads from a metagenomic dataset using an overlap graph. Note that OGRE is not fully developed and optimized
software. However the code is available and usable for read clustering.

---------------
 Installation
---------------

1. Download OGRE: https://github.com/Marleen1/OGRE

2. Install the following dependencies:
- Python 3 with packages numpy, scikit-learn and scipy
- Snakemake (https://snakemake.readthedocs.io/en/stable/)
- Minimap2 (https://github.com/lh3/minimap2)

--------------------
 Setting the inputs
--------------------

The header of the file "Snakefile" contains the required input information. The following need to be specified in the file header:
- folder: the location of the read data.
- fqfile: the filename of the fastq file containing the reads that need to be clustered.
- limits: the desired maximum cluster size, given as a string. For example a maximum cluster size of 5000 should be written as '5000'. A list of cluster sizes can be
provided as well, OGRE will then output multiple clusterings for various maximum cluster sizes. For example: limits = ['5000','10000','inf'].
- numsplitlines: OGRE starts by splitting up the original fastq file into smaller files. This is the number of lines in each small fastq file. If you have a large amount
of RAM, you can set this value high. On a system with 125 GB RAM we used numsplitlines = 8000000.
- train_data: the training dataset. A training dataset is provided with OGRE. In order to use this dataset, choose train_data = 'traindata_CAMI23toy23_80000.txt'
- lines_per_file: OGRE creates an overlap file, which represents the overlap graph. This file is later split into smaller files and processed. If you have a large amount
of RAM, you can choose a high value here. we chose lines_per_file = 2500000 on a computer with 125 GB RAM.

--------------
 Running OGRE
--------------

OGRE is run by typing

snakemake -j xx

in the command line, where xx is the number of threads it can use.
