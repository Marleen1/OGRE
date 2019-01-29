# OGRE
OGRE is an overlap graph-based read assembly tool.

All code for the development of OGRE is provided here.

OGRE assumes the following input files:
- a fastq file of interleaved paired-end reads
- a training data file

# Setting up
The heading of Snakefile contains all its inputs.

Choose one of the CAMI datasets that you wish to cluster. This determines which of the traindata files to use: pick the one where the data to be clustered is not included. For example, for clustering the CAMI_low dataset (or CAMI1), use traindata_CAMI23_toy12_80000.txt as the training dataset.

Please provide the following parameters in the header of Snakefile:
- numsplits determines the number of sub fastq files the original fastq file is split into for the overlap graph construction step (using Minimap2). Using more splits requires a longer processing time, but limits the required amount of free disk space;
- folder is the folder containing the data (so the fastq file);
- fqfile is the name of the fastq file, without the ".fq" extension;
- numlines is the number of lines of the fastq file;
- lines_per_file is the number of lines in each split file for the prediction step;
- subset is the data subset (only applies to CAMI_medium, CAMI_high, toy_medium and toy_high);
- train_data is the training dataset;
- val_data is the dataset that will be clustered (denoted as CAMI1, CAMI2, CAMI3, toy2 or toy3);
- run_ID is the ID that will be used in all output files.

