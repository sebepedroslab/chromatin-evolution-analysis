#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h
#$ -l virtual_free=50G,h_rt=43200
#$ -o tmp/
#$ -e tmp/

#Defining arguments
if [ -z "$3" ]
	then echo -e "\nERROR! You need to specify arguments!

	- 1. input fasta
	- 2. Pfam database (hmm repository)
	- 3. num cpu
"
    exit
fi

# function
function do_pfamscan {

    # input
    fas=$1
    pfam_db=$2
    n_cpu=$3

    # pfamscan
    pfam_scan.pl -fasta ${fas}  -dir ${pfam_db} -e_seq 0.01 -e_dom 0.01 > ${fas%%.fasta}.pfamscan.csv


    # list of proteins
    bioawk -c fastx '{ print $1 }' ${fas} > ${fas}_TMP_list

    # table with pfam architectures
    while read seqid ; do
    grep -w "$seqid" ${fas%%.fasta}.pfamscan.csv \
    | tr -s  ' ' '\t' | cut -f 7 | xargs
    done < ${fas}_TMP_list > ${fas}_TMP_architectures
    paste ${fas}_TMP_list ${fas}_TMP_architectures > ${fas%%.fasta}.pfamscan_archs.csv

    # clean
    rm ${fas}_TMP_architectures ${fas}_TMP_list

}

do_pfamscan $1 $2 $3

