

prefix=inf
target=inf
dataset=294_220 
clusters_json=/export/scratch2/vincent/project/OGRE/data/$dataset/"$dataset"_sr_k21_w11_s60_m60_n2_r0_A4_B2_max$prefix"_final_clusters_grouped.json"


python  evaluate_clusters.py $dataset.$target.from.$prefix  $dataset.$target.from.$prefix.clustersize.json  /export/scratch2/vincent/project/OGRE/3.mash/$dataset

echo "all finished..."

