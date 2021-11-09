#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=60G,h_rt=2592000
#$ -o tmp/
#$ -e tmp/

# input
i=$1
c=$2
o=$3

#mafft --genafpair --thread $c --reorder --maxiterate 10000 $i > ${i%%.fasta}.l.fasta
#clipkit ${i%%.fasta}.l.fasta -m kpic-gappy -o ${i%%.fasta}.lt.fasta -g 0.7
/users/asebe/xgraubove/Programes/iqtree-2.1.0-Linux/bin/iqtree2 -s ${i%%.fasta}.lt.fasta -m LG+G4 -nt AUTO -ntmax $c -pre ${o} -nm 2000 -cptime 1800 -pers 0.1 -bb 1000 -te euk.Histone.iqt.treefile

# now check if there are outlier sequences in the tree
python /users/asebe/xgraubove/metazoan-tf-phylo/scripts/run_treeshrink.py -c -t ${o}.treefile -m per-gene -q 0.01 -s 10,1 -f ${o}.treeshrink.firstpass.txt

