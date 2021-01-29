
#1.get fastq for each cluster
less cluster2reads.json |sed 's/"/\n/g'|grep '^[0-9]' >cluster.id

split_file cluster.id 50 cluster.id_split


perl -e 'for my$i(10..17){print "python get_fq4cluster_subfile.py all.fq cluster2reads.json cluster.id_split_$i true $i\n";}' >run.sh

#submit jobs
#qsub.pl -l h_rt=03:00:00,vf=6G -m 80 -s 60 run.sh &

cat run.sh|xargs -i -P 16 bash -c "{}" &

