# copy files first: *.py *.sh *.R  read_species_dict.json

prefix=inf
target=inf
dataset=294_220
clusters_json=/export/scratch2/vincent/project/OGRE/2.fq4cluster/294_220/inf/cluster2reads.json 


for i in /export/scratch2/vincent/project/OGRE/2.fq4cluster/$dataset/$prefix/fastq_*;do export ii=$i;ls $i|sed 's/\///'|perl -ne 'chomp;print "$ENV{ii}/$_/$_.fq\n";';done >fq.list

mash sketch -p 20 -l fq.list -o sketch.out
echo "skethc finished..."
#remove p value=1 because the outfile is too big
mash dist -p 20 sketch.out.msh sketch.out.msh |perl -ne 'my@a=split; next if $a[0] eq $a[1];next if $a[3]==1;$a[0]=~s/\S+\///;$a[0]=~s/\.fq$//;$a[1]=~s/\S+\///;$a[1]=~s/\.fq$//;my$min=$a[0];my$max=$a[1];if($a[0]>$a[1]){$min=$a[1];$max=$a[0];} print "$min\t$max\t$a[2]\t$a[3]\n";' > dist.out

cat dist.out |perl -ne 'BEGIN{%h;}my@a=split; my$min=$a[0];my$max=$a[1];next if exists $h{"$min:$max"};$h{"$min:$max"}=1; print "$min\t$max\t$a[2]\t$a[3]\n";'  >distances.tab

echo "dist finished..."
rm sketch.out.msh
#rm dist.out

cut -f 4 distances.tab >p.txt

Rscript padj.R
echo "p value adjusting finished..."

paste distances.tab padj.txt >distances.padj.tab

rm -f distances.tab p.txt padj.txt

##

cat distances.padj.tab|awk '$5<=0.01' |sort -gk3 >distances.padj.tab.filter

python ./merge_clusters_from_mash.v2.py $clusters_json distances.padj.tab.filter $target $dataset.$target.from.$prefix.clusters_grouped.json

# evaluate

python get_cluster_size_json.py  $dataset.$target.from.$prefix.clusters_grouped.json  $dataset.$target.from.$prefix.clustersize.json


sh evaluate_clusters.sh



cat hist_species_per_cluster_*.txt


echo "all finished..."

