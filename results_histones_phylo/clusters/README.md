# Histone clusters

Steps:

1. Retrieve assembly-level coordinates for each histone, and their closest downstream and upstream genes:

```bash
s01_txvalidation_assembly_2020-12-25.sh
```

2. Annotate histone classes, and find clusters:

```bash
Rscript s02_annotate_histone_clusters.R
```

3. Run [EvolClust](https://github.com/Gabaldonlab/EvolClust/wiki/2.--Running-Evolclust) to identify clusters of homology, to ascertain whether histones can be found in them:

```bash
### Install Evolclust:
# git clone git@github.com:Gabaldonlab/EvolClust.git
# conda create -n evolclust python=2.7 numpy
# conda activate evolclust

# obtain evolclust-friendly gene IDs for the species of interest:
Rscript s03_evolclust_prepare_data.R

# run orthofinder with low inflation, to obtain big OGs:
# (only for species of interest that have H2B-H2A or H3-H4 clusters)
cd /users/asebe/xgraubove/histonome-ops/results-orthofinder-small
qsub -N ofhist -pe smp 12 qsub_orthofinder.sh input_fasta/ 1.2 12

# format OF output to fit evolclust:
cut -f 2- -d ':' /users/asebe/xgraubove/histonome-ops/results-orthofinder-small/input_fasta/OrthoFinder/Results_Mar14/Orthogroups/Orthogroups.txt | sed "s/^ //"| tr ' ' '\t' > data-evolclust/homology.mcl
# use proper gene IDs:
awk 'NR==FNR { l[$1]=$2;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' <(cat data-evolclust/*.gene_dict.csv) data-evolclust/homology.mcl | tr ' ' '\t' > data-evolclust/homology-evolclust.mcl
tr '\t' '\n' < data-evolclust/homology-evolclust.mcl > data-evolclust/homology-evolclust.list.txt

# run evolclust:
#conda activate evolclust
#change dit to evolclust local folder????
python evolclust.py -i /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_histones_phylo/clusters/data-evolclust/homology-evolclust.mcl -l /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_histones_phylo/clusters/data-evolclust/homology-evolclust.list.txt -d results-evolclust --local
# QUITE BIG, DO IT IN THE SAME FOLDER AS ORTHOFINDER

# post-processing in R:
# - load histonome-to-assembly gene mappings
# - load assembly-to-evolclust mappings
# - load evolclust clusters
# - check whether we can find histones within any of the clusters, and the phylogenetic extent of said clusters
```
