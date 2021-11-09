#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=10G,h_rt=2592000
#$ -o tmp/
#$ -e tmp/

# function
function do_pfamscan {

    # input
    fas=$1
    pfam_db=$2
    n_cpu=$3

    # pfamscan
    pfam_scan.pl  -cpu ${n_cpu}  -fasta ${fas}  -dir ${pfam_db} > ${fas%%.fasta}.pfamscan.csv

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

