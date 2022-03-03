#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h,long-sl7,short-sl7
#$ -l virtual_free=30G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/

bs=$1
ou=$2
on=$3

# output folder
mkdir -p ${ou}

echo "output folder: ${ou}"
ls ${ou}
echo "input files: ${bs}"

# # index output
samtools merge -f ${ou}/${on}.bam $(echo ${bs} | tr ',' ' ') -@ 10
samtools index ${ou}/${on}.bam

# only primary alignments
#samtools view -bS -F 256 ${ou}/${on}.bam > ${ou}/${on}.best.bam
#samtools index ${ou}/${on}.best.bam

# # bedgraph coverage
bedtools genomecov -ibam ${ou}/${on}.bam -bga > ${ou}/${on}.bdg
#bedtools genomecov -ibam ${ou}/${on}.best.bam -bga > ${ou}/${on}.best.bdg

# counts
samtools idxstats ${ou}/${on}.bam > ${ou}/${on}.counts.tsv
#samtools idxstats ${ou}/${on}.best.bam > ${ou}/${on}.best.counts.tsv
