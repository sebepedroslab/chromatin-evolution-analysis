# Data folders

* `data`: original data and ubiquitin searches. These files will be parsed from the corresponding Excel files.
* `data-paired-searches`: data from paired searches of acetylation with mono/di/trimethylation (three files per species). These files are parsed via `sqlite3` (see `scripts` folder).

# Retrieve predicted peptides from prokaryotes

Get sequences:

```bash
conda activate base
# select species
i="Methanobrevibacter cuticularis"
i="Nitrososphaera viennensis"

# get sequences
mkdir -p ~/histonome-ops/data/Proteomes_proteomics_identification-2021/fasta_named/
mkdir -p ~/histonome-ops/data/Proteomes_proteomics_identification-2021/fasta/
j=$(echo $i | sed "s/ /_/" )
bioawk -c fastx '{ print $1,$2,$4 }' ~/histonome-ops/data/sequences/seq_Archaea.fasta | grep -w "${i}" | awk '{ print ">'${j}'_"$1"\n"$2 }' >  ~/histonome-ops/data/Proteomes_proteomics_identification-2021/fasta/${j}.fasta

# tag histones manually
bioawk -c fastx '{ print $1 }' ..//results_histones_phylo/archaeal_Ntails/arc.Histone.seqs.fasta| fgrep -f - ~/histonome-ops/data/Proteomes_proteomics_identification-2021/fasta/${j}.fasta
nano  ~/histonome-ops/data/Proteomes_proteomics_identification-2021/fasta/${j}.fasta # add "|Histone" tag manually
```
