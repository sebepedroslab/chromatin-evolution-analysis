#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=100G,h_rt=720:00:00
#$ -o tmp/
#$ -e tmp/

# input
i=$1 # NON-ALIGNED FASTA
c=$2 # NUM CPUS

function do_alitrimphy {

	local i=$1
	local c=$2
	local il=$3
	local ilt=$4
	local tre=$5

	# align & trim
	if [ -s $il ] ; then
	echo "$il already exists, go to next step"
	else
	mafft --globalpair --thread $c --reorder --maxiterate 10000 $i > $il
	fi
	if [ -s $ilt ] ; then
        echo "$ilt already exists, go to next step"
	else
	#trimal -in ${i%%.fasta}.l.fasta -out ${i%%.fasta}.lt.fasta -gappyout
	clipkit $il -m kpic-gappy -o $ilt -g 0.7
	fi

	# iqtree
	/users/asebe/xgraubove/Programes/iqtree-2.1.0-Linux/bin/iqtree2 \
		-s $ilt \
		-m TEST \
		-mset LG,WAG,JTT \
		-nt AUTO \
		-ntmax $c \
		-bb 1000 \
		-pre $tre \
		-nm 10000 \
		-nstop 200 \
		-cptime 1800

}

# do alignments and trees
do_alitrimphy $i $c ${i%%.fasta}.l.fasta ${i%%.fasta}.lt.fasta ${i%%.fasta}.iqtree



