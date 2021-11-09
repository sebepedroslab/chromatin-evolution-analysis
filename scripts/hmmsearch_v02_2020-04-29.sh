#!/bin/bash

# Defining arguments
if [ -z "$6" ]
	then echo -e "
ERROR! You need to specify arguments!
- 1 hmm file
- 2 output prefix (can include folder path)
- 3 database to search (FASTA, must include blastdb)
- 4 PFAM-A db folder (pfam_scan.pl). If NA, pfamscan is skipped
- 5 evalue (GA or 1e-5, etc)
- 6 num cpus (hmmsearch and pfam_scan.pl)
"
	exit
fi

f="$1"
jobid="$2"
answerDB="$3"
answerPfam="$4"
ev="$5"
n_cpu="$6"
if [ ! ${answerPfam} = "NA" ] ; then pfam_db=${answerPfam} ; fi

#Look for hmmsearch hits, extract hits, computes architectures, computes taxonomic report

echo -e "\n>>> Looking for $f with hmmsearch in $answerDB <<<\n"

if [ $ev == GA ] ; then

	hmmsearch \
		-A ${jobid}_TMP_align.stockholm \
		--cut_ga \
		--cpu ${n_cpu} \
		${f} \
		${answerDB} \
		> ${jobid}.hmmerout

else

	hmmsearch \
		-A ${jobid}_TMP_align.stockholm \
		-E ${ev} \
                --cpu ${n_cpu} \
		--incE ${ev} \
		--incdomE ${ev} \
		${f} \
		${answerDB} \
		> ${jobid}.hmmerout

fi

# extract domains
esl-reformat fasta ${jobid}_TMP_align.stockholm \
	| sed "s/\//_/" \
	| bioawk -c fastx ' { print ">"$1"\n"$2 }' \
	> ${jobid}_domains.fasta

# sequences list
esl-reformat fasta ${jobid}_TMP_align.stockholm \
	| bioawk -c fastx '{ print $1 }' \
	| cut -f1 -d '/' \
	| sort -u > ${jobid}_TMP_hits.list

# extract whole sequences
blastdbcmd \
	-db $answerDB \
	-dbtype prot \
	-entry_batch ${jobid}_TMP_hits.list \
	-outfmt "%f" \
	> ${jobid}_seqs.fasta

# report
numhits=$(grep -c ">" ${jobid}_seqs.fasta)
numdoms=$(grep -c ">" ${jobid}_domains.fasta)
echo -e "$numhits seqs extracted in general fasta"
echo -e "$numdoms seqs extracted in domains fasta"
if [ $numhits -eq 0 ] ; then echo "no hits! exit" ; fi
if [ $numdoms -eq 0 ] ; then echo "no domains! exit" ; fi


# Pfamsearch (opcional)
if [ ! -z ${pfam_db} ] ;
then echo -e "Computing domain architectures with Pfamscan..."
pfam_scan.pl \
	-cpu ${n_cpu} \
	-fasta ${jobid}_seqs.fasta \
	-dir ${pfam_db} \
	> ${jobid}_seqs.pfamscan
	while read seqid ; do
		grep "$seqid" ${jobid}_seqs.pfamscan \
		| tr -s  ' ' '\t' \
		| cut -f 7 \
		| xargs
	done < ${jobid}_TMP_hits.list > ${jobid}_TMP_PFAMTABLE_RESULTS
	paste ${jobid}_TMP_hits.list ${jobid}_TMP_PFAMTABLE_RESULTS > ${jobid}_seqs.domain_architectures.csv
else echo -e "Skiping Pfamscan search..."
fi

#Ends...
echo -e "\n\nRemoving temporary files..."
rm -f ${jobid}_TMP_*

echo -e "\n\n>> DONE! <<\n\n"

exit

