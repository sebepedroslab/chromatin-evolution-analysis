# Classify canonical histones

Recipe

## Classification of histones

1. List of reference histones from human (Ensembl Release 100): `human_histones.dict_classes.csv`.

2. Pairwise alignments to human reference:

```bash
# reference dbs
diamond makedb --in histone_db_ref.fasta -d histone_db_ref.fasta
diamond makedb --in human_histones.dict_classes.fasta -d human_histones.dict_classes.fasta

# euks
diamond blastp -d histone_db_ref.fasta              -q euk.Histone.domains.fasta -o euk.Histone.to_histdb.csv --more-sensitive --max-target-seqs 100 --quiet && gzip euk.Histone.to_histdb.csv
diamond blastp -d human_histones.dict_classes.fasta -q euk.Histone.domains.fasta -o euk.Histone.to_human.csv  --more-sensitive --max-target-seqs 100 --quiet && gzip euk.Histone.to_human.csv

# add virus and archaea too
cat arc.Histone.domains_wp.fasta vir.Histone.domains_wp.fasta euk.Histone.domains.fasta > all.Histone.domains.fasta
diamond blastp -d histone_db_ref.fasta              -q all.Histone.domains.fasta -o all.Histone.to_histdb.csv --more-sensitive --max-target-seqs 100 --quiet && gzip all.Histone.to_histdb.csv
diamond blastp -d human_histones.dict_classes.fasta -q all.Histone.domains.fasta -o all.Histone.to_human.csv  --more-sensitive --max-target-seqs 100 --quiet && gzip all.Histone.to_human.csv
```

3. Pairwise alignments to self:

```bash
# eukaryota
diamond makedb --in euk.Histone.domains.fasta -d euk.Histone.domains.fasta
diamond blastp -d euk.Histone.domains.fasta -q euk.Histone.domains.fasta -o euk.Histone.to_self.csv --more-sensitive --max-target-seqs 100 --quiet && gzip euk.Histone.to_self.csv

# euk+vir+archaea
diamond makedb --in all.Histone.domains.fasta -d all.Histone.domains.fasta
diamond blastp -d all.Histone.domains.fasta -q all.Histone.domains.fasta -o all.Histone.to_self.csv --more-sensitive --max-target-seqs 100 --quiet && gzip all.Histone.to_self.csv
```

4. Alignments to histones used for proteomics analyses:

```bash
diamond makedb --in all_proteomes_for_proteomics_2021-09-08.fasta -d ~/histonome-ops/data/all_proteomes_for_proteomics_2021-09-08.fasta
diamond blastp -d all_proteomes_for_proteomics_2021-09-08.fasta -q all.Histone.domains.fasta -o all.Histone.to_proteomics.csv --more-sensitive --max-target-seqs 100 --quiet && gzip all.Histone.to_proteomics.csv
```

5. Create histone networks:

```bash
# networks from eukaryotic and all alignments to histdb
python s01_project_network_alignments-naive.py euk.Histone.to_self.csv.gz euk.Histone.to_histdb.csv.gz alignments
python s01_project_network_alignments-naive.py all.Histone.to_self.csv.gz all.Histone.to_histdb.csv.gz alignments

# additional networks to human (same results, but TFs look better)
python s01_project_network_alignments-naive.py euk.Histone.to_self.csv.gz euk.Histone.to_human.csv.gz alignments
python s01_project_network_alignments-naive.py all.Histone.to_self.csv.gz all.Histone.to_human.csv.gz alignments

# explore annotations a bit with this R script:
Rscript s02_classification_from_human_ref.R
```

5. Create composite of connected components:

```bash
# to obtain bitmaps
for i in all.Histone.to_histdb.Spring-Naive-cla.*.ne.svg ; do inkscape --with-gui --select=LineCollection_1 --verb SelectionCreateBitmap $i & \ ; done
killall inkscape

# To create composites from PDF:
# pdfjam --nup 4x3 euk.Histone.to_histdb.Spring-Naive-cla.pdf --outfile ../euk.Histone.to_histdb.Spring-Naive-composite.pdf
# pdfjam --nup 5x5 all.Histone.to_histdb.Spring-Naive-cla.pdf --outfile ../all.Histone.to_histdb.Spring-Naive-composite.pdf
```

## Targeted phylogenies

1. Extract sequences in each connected component of the networks above, and align them:

```bash
# eukaryota only groups
for i in alignments/euk.Histone.to_histdb/CC*.txt ;    do xargs faidx -d ' ' euk.Histone.domains.fasta < $i >  ${i%%.txt}.fasta ; mafft --genafpair --thread 1 --reorder --maxiterate 1000 ${i%%.txt}.fasta > ${i%%.txt}.l.fasta ; done
for i in alignments/euk.Histone.to_histdb/CC*.l.fasta; do iqtree2 -s ${i} -m TEST -mset LG -nt AUTO -ntmax 8 -fast -pre ${i%%.l.fasta}.iqt ; done

# include viruses and archaea
for i in alignments/all.Histone.to_histdb/CC*_H*.txt ; do 
xargs faidx -d ' ' all.Histone.domains.fasta < $i >  ${i%%.txt}.fasta
mafft --genafpair --thread 4 --reorder --maxiterate 1000 ${i%%.txt}.fasta > ${i%%.txt}.l.fasta
iqtree2 -s ${i%%.txt}.l.fasta -m TEST -mset LG -nt AUTO -ntmax 8 -fast -pre ${i%%.txt}.iqt
done

# plot histone trees
Rscript s03_paint_phylogenies.R

# investigate relative distances between canonical, non canonical etc. histones in each tree
Rscript s04_diversification_euk_histones.R
```

2. Concatenate all main components that contain viral sequences, and align them:

```bash
# which components contain viral sequences?
# RUN: grep "^vir" all.Histone.to_histdb.Spring-Naive.csv | cut -f2 | sort | uniq -c| sort -n
      2 37 # nudivirus H3-H4 (single seq)
      3 38 # TF?? exclude
      5 33 # H2B
      6 41 # H4-like, giant virus
      7 30 # H3
     12 39 # H2A-like, giant virus only
     12 40 # H3-like, giant virus only
     13 34 # H2A
     15 31 # H4
```

```bash
# select key components:
for i in 37 33 41 30 39 40 34 31 ; do 
awk '$2 == "'$i'" { print $1 }' all.Histone.to_histdb.Spring-Naive.csv | xargs faidx -d ' ' all.Histone.domains.fasta
done > alignments/vir.concatenated/allhistone_dom.fasta

# cdhit? probably not necessary, doesn't significantly solve the problem
# cd-hit -i alignments/vir.concatenated/allhistone_dom.fasta -o alignments/vir.concatenated/allhistone_dom.c1.fasta -c 1 -d 0

# align:
mafft --globalpair --thread 4 --reorder --maxiterate 1000 alignments/vir.concatenated/allhistone_dom.fasta > alignments/vir.concatenated/allhistone_dom.g.fasta
clipkit alignments/vir.concatenated/allhistone_dom.g.fasta -m kpic-gappy -o alignments/vir.concatenated/allhistone_dom.gt.fasta -g 0.7
# tree
# iqtree2 -s alignments/vir.concatenated/allhistone_dom.g.fasta -m TEST -nt AUTO -ntmax 8 -pre alignments/vir.concatenated/allhistone_dom.iqtfb -nm 10000
iqtree2 -s alignments/vir.concatenated/allhistone_dom.gt.fasta -m TEST -nt AUTO -ntmax 8 -pre alignments/vir.concatenated/allhistone_dom.gt.iqtfb -nm 10000 -bb 1000
iqtree2 -s alignments/Eri -m TEST -nt AUTO -ntmax 8 -pre alignments/vir.concatenated/allhistone_dom.gt.iqtfb -nm 10000 -bb 1000



#### ALTERNATIVE:
# add hits from Erives 2017 to have good representation of Marseillevirus sequences and an archaeal outgroup
# manually add hits from Erives 2017 phylogenies:
# first, retrieve blast hits:
diamond blastp -d ~/histonome-ops/data/sequences/seq_Virus.fasta -q alignments/Erives_2017.viral_histones.fasta -o alignments/Erives_2017.viral_histones.to_NCBI.diamond.csv --more-sensitive --max-target-seqs 100 --quiet
# get hits and their coordinates, and remove redundancy in hits coordinates
awk 'BEGIN {OFS="\t" } { print $2,$9,$10 }' alignments/Erives_2017.viral_histones.to_NCBI.diamond.csv > alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits.txt
bedtools merge -d -10 -i <(sort -k 1,1 -k 2,3n alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits.txt ) > alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.txt
esl-sfetch -C -f ~/histonome-ops/data/sequences/seq_Virus.fasta <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.txt ) > alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.fasta
sed -i "s/>/>vir_/" alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.fasta
sed -i "s/ /|Erives2017 /" alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.fasta

# manually add a small archaeal outgroup based on Erives 2017:
diamond blastp -d ~/histonome-ops/data/sequences/seq_Archaea.fasta -q alignments/Erives_2017.archaeal_histones.fasta -o alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.csv --more-sensitive --max-target-seqs 100 --quiet
# get hits and their coordinates, and remove redundancy in hits coordinates
awk 'BEGIN {OFS="\t" } { if( $8>50) { print $2,$9,$10 }  }' alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.csv > alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits.txt
bedtools merge -d -10 -i <(sort -k 1,1 -k 2,3n alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits.txt ) > alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.txt
esl-sfetch -C -f ~/histonome-ops/data/sequences/seq_Archaea.fasta <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.txt ) > alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.fasta
sed -i "s/>/>arc_/" alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.fasta
sed -i "s/ /|Erives2017 /" alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.fasta
cd-hit -i alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.fasta -o alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.c95.fasta -c 0.95 -d 0

# add viruses and archaea into original dataset:
cat alignments/vir.concatenated/allhistone_dom.fasta > alignments/vir.concatenated/cat_histones.fasta
cat alignments/Erives_2017.viral_histones.to_NCBI.diamond.hits_nr.fasta >> alignments/vir.concatenated/cat_histones.fasta
cat alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.c95.fasta >> alignments/vir.concatenated/cat_histones.fasta

# run phylogeny
qsub -N virhis -pe smp 8 qsub_alignment-alrt.sh alignments/vir.concatenated/cat_histones.fasta 8


#### CRAZY IDEA
# add Ribosomal_S6e homologs from eukaryotes and archaea
cp ~/histonome-ops/results-phylogenies/searches_HistRib_Abr21/euk.Ribosomal_S6e.domains.fasta ribo.euk.domains.fasta
cd-hit -i ~/histonome-ops/results-phylogenies/searches_HistRib_Abr21/arc.Ribosomal_S6e.domains.fasta -o ribo.arc.domains.fasta -c 0.95 -d 0
mkdir -p alignments/ribo.concatenated
cat alignments/vir.concatenated/cat_histones.fasta ribo.euk.domains.fasta ribo.arc.domains.fasta > alignments/ribo.concatenated/cat_ribohist.fasta
qsub -N virribohis -pe smp 8 qsub_alignment-alrt.sh alignments/ribo.concatenated/cat_ribohist.fasta 8
```

## Old

Concatenate all main components and then align them:

* include archaea and viruses
* most euk-like archaeal histones belong to CC0
* histone dimers from Euryarchaeota seem to belong to CC4

```bash
# concatenate alignments, store CC, remove dashes
for i in alignments/all.Histone.to_histdb/*.l.fasta ; do sed "s/\([0-9]\)$/\1|$(basename $i | sed "s/.l.fasta//" | sed "s/_$//" )/" $i ; done | bioawk -c fastx '{ gsub(/-/,"",$2) ; print ">"$1"\n"$2 }' > alignments/all.concatenated/allhistone_dom.fasta

# cdhit -- REALLY NECESSARY!
cd-hit -i alignments/all.concatenated/allhistone_dom.fasta -o alignments/all.concatenated/allhistone_dom.c99.fasta -c 0.99 -d 0

# align
mafft --globalpair --thread 8 --reorder --maxiterate 1000 alignments/all.concatenated/allhistone_dom.fasta > alignments/all.concatenated/allhistone_dom.g.fasta
# tree
iqtree2 -s alignments/all.concatenated/allhistone_dom.g.fasta -m TEST -nt AUTO -ntmax 8 -pre alignments/all.concatenated/allhistone_dom.iqtfb -nm 10000

# alternative: only main components (not a great idea, some important archaea are missing)
# mkdir alignments_concatenated
# # concatenate
# cat alignments/all.Histone.to_histdb/CC0_.fasta alignments/all.Histone.to_histdb/CC51_H2A.fasta alignments/all.Histone.to_histdb/CC50_H2B.fasta alignments/all.Histone.to_histdb/CC46_H3.fasta alignments/all.Histone.to_histdb/CC47_H4.fasta  > alignments_concatenated/archeuk_histones.fasta
# # cdhit 
# cd-hit -i alignments_concatenated/archeuk_histones.fasta -o alignments_concatenated/archeuk_histones.c099.fasta -c 0.99 -d 0
# # align and tree
# mafft --genafpair --thread 1 --reorder --maxiterate 1000 alignments_concatenated/archeuk_histones.c099.fasta > alignments_concatenated/archeuk_histones.c099.l.fasta
# iqtree2 -s alignments_concatenated/archeuk_histones.c099.l.fasta -m TEST -mset LG -nt AUTO -ntmax 8 -fast -pre alignments_concatenated/archeuk_histones.c099.iqt
```

3. Investigate diversification in the main canonical components:

```bash
Rscript
```
