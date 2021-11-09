#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m abe
#$ -q long-sl7
#$ -l virtual_free=10G,h_rt=2592000
#$ -o tmp/
#$ -e tmp/

makeblastdb -dbtype prot -parse_seqids -in data/databases_seqs/seq_Bacteria.fasta
