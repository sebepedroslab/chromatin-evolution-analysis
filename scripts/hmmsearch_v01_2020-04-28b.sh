#!/bin/bash

# Defining arguments
if [ -z "$7" ]
then echo -e "
ERROR! You need to specify arguments!
- 1 hmm file
- 2 output prefix
- 3 whether you want pfamscan or not (y/n)
- 4 database to search (FASTA, must include blastdb)
- 5 sps list (només extrau aquestes)
- 6 num cpus (hmmsearch and pfam_scan)
- 7 PFAM-A db folder (pfam_scan)
"
exit
fi

f="$1"
jobid="$2"
answerPfam="$3"
answerDB="$4"
sl="$5"
n_cpu="$6"
pfam_db="$7"

#Look for hmmsearch hits, extract hits, computes architectures, computes taxonomic report

echo -e "\n>>> Looking for $f with hmmsearch in $answerDB <<<\n"

hmmsearch \
	--cut_ga \
	--tblout ${jobid}.tblout \
	--domtblout ${jobid}.domtblout \
	--pfamtblout ${jobid}.pfamtblout \
	--cpu ${n_cpu} \
	${f} \
	${answerDB} \
	> ${jobid}.hmmerout


# prepare species list for taxonomic reports?
sed -e "s/$/_/" \
	-e "s/^/[subseq from] /" ${sl} > ${jobid}_TMP_pseudolist

sed -e "s/$/_/" \
	-e "s/^/>/" ${sl} > ${jobid}_TMP_pseudolist2

fgrep -f ${jobid}_TMP_pseudolist ${jobid}_TMP_align.stockholm \
	| tr -s ' ' '\t' \
	| tr -s '/' '\t' \
	| cut -f 2 \
	| sort -u \
	> ${jobid}_TMP_log

esl-reformat \
	fasta \
	${jobid}_TMP_align.stockholm \
	| sed "s/\[subseq from\] .*//" \
	| tr "/" " " \
	| bioawk -c fastx '{ print ">"$1,$4"	"$2 }' \
	| fgrep -f ${jobid}_TMP_pseudolist2 \
	| awk '{ print $1,$2"\n"$3 }' \
	> ${jobid}_domains.fasta

grep ">" ${jobid}_domains.fasta \
	| sed "s/>//" \
	| cut -f1 -d ' ' \
	| sort \
	| uniq -c \
	| awk '$1 > 1 {print $2,$1}' > ${jobid}_domains.multiplereport

numseqs=$(sort -u ${jobid}_TMP_log \
		| sed '/^$/d' \
		| wc -l \
		| sed -e 's/^[ \t]*//')

echo -e "$numseqs seqs found by hmmscan"

blastdbcmd \
	-db $answerDB \
	-dbtype prot \
	-entry_batch ${jobid}_TMP_log \
	-outfmt "%f" \
	> ${jobid}_seqs.fasta

numhits=$(grep -c ">" ${jobid}_seqs.fasta)
numdoms=$(grep -c ">" ${jobid}_domains.fasta)

echo -e "$numhits seqs extracted in general fasta"
echo -e "$numdoms seqs extracted in domains fasta"


echo -e "Taxonomic report..."

while read si ; do
	echo "$si $(grep -c ">${si}_" ${jobid}_seqs.fasta)"
done < ${sl} > ${jobid}_seqs.taxonreport

echo -e "\nMultiple domain report...
($(wc -l ${jobid}_domains.multiplereport))"
head -n 20 ${jobid}_domains.multiplereport
echo "..."

# Pfamsearch (opcional)

if [ $answerPfam = y ] ;
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
		done < ${jobid}_TMP_log > ${jobid}_TMP_PFAMTABLE_RESULTS
		paste ${jobid}_TMP_log ${jobid}_TMP_PFAMTABLE_RESULTS > ${jobid}_seqs.domain_architectures.csv
	else echo -e "Skiping Pfamscan search..."
fi

#Ends...
echo -e "\n\nRemoving temporary files..."
# rm -f ${jobid}_TMP_*

echo -e "\n\n>> DONE! <<\n\n"

exit

