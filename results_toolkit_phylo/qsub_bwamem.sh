#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h,long-sl7,short-sl7
#$ -l virtual_free=10G,h_rt=6:00:00
#$ -o tmp/
#$ -e tmp/

rl=$1
r1=$1
r2=$2
ix=$3
th=$4
ou=$5

# create output folder
mkdir -p ${ou}

# run bwa mem
bwa mem \
	-M -t ${th} \
	${ix} \
	${r1} \
	${r2} \
	| bamsormadup inputformat=sam threads=${th} tmpfile=${ou}/$(basename ${r1%%_1.fastq.gz}).tmp.$(date +%s) SO=coordinate \
	> ${ou}/$(basename ${r1%%_1.fastq.gz}).bwa.bam

# index output	
samtools index ${ou}/$(basename ${r1%%_1.fastq.gz}).bwa.bam

# bedgraph coverage
bedtools genomecov -ibam ${ou}/$(basename ${r1%%_1.fastq.gz}).bwa.bam -bga > ${ou}/$(basename ${r1%%_1.fastq.gz}).bwa.bdg
