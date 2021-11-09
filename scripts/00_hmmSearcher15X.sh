#!/bin/bash

# Runs in a folder with .hmm profiles

rm -f *_TMP_*

#Defining arguments
if [ -z "$6" ]
	then echo -e "\nERROR! You need to specify arguments!
- 1 hmm file
- 2 job id
- 3 whether you want pfamscan or not (y/n)
- 4 database to search:
- 5 sps list (només extrau aquestes)
- 6 evalue
    GA = use Gathering Threshold
    1e-5  = use 1e-5
	"
	exit
fi

f="$1"
jobid="$2"
answerPfam="$3"
answerDB="$4"
sl="$5"
ev="$6"


#Look for hmmsearch hits, extract hits, computes architectures, computes taxonomic report

echo -e "\n>>> Looking for $f with hmmsearch in $answerDB <<<\n"

if [ $ev == GA ] ; then

	hmmsearch \
		-A ${f%%.*}_${jobid}_TMP_align.stockholm \
		--cut_ga \
		$f \
		$answerDB \
		> ${f%%.*}_${jobid}.hmmerout

else

	hmmsearch \
		-A ${f%%.*}_${jobid}_TMP_align.stockholm \
		-E $ev \
		--incE $ev \
		--incdomE $ev \
		$f \
		$answerDB \
		> ${f%%.*}_${jobid}.hmmerout


fi

sed -e "s/$/_/" \
	-e "s/^/[subseq from] /" $sl > _TMP_pseudolist

sed -e "s/$/_/" \
	-e "s/^/>/" $sl > _TMP_pseudolist2


fgrep -f _TMP_pseudolist ${f%%.*}_${jobid}_TMP_align.stockholm \
	| tr -s ' ' '\t' \
	| tr -s '/' '\t' \
	| cut -f 2 \
	| sort -u \
	> ${f%%.*}_${jobid}_TMP_Nyeidlist.log
	

esl-reformat \
	fasta \
	${f%%.*}_${jobid}_TMP_align.stockholm \
	| sed "s/\[subseq from\] .*//" \
	| tr "/" " " \
	| bioawk -c fastx '{ print ">"$1,$4"	"$2 }' \
	| fgrep -f _TMP_pseudolist2 \
	| awk '{ print $1,$2"\n"$3 }' \
	> ${f%%.*}_${jobid}_domains.fasta

grep ">" ${f%%.*}_${jobid}_domains.fasta \
	| sed "s/>//" \
	| cut -f1 -d ' ' \
	| sort \
	| uniq -c \
	| awk '$1 > 1 {print $2,$1}' > ${f%%.*}_${jobid}_domains.multiplereport

numseqs=$(sort -u ${f%%.*}_${jobid}_TMP_Nyeidlist.log \
		| sed '/^$/d' \
		| wc -l \
		| sed -e 's/^[ \t]*//')

echo -e "$numseqs seqs found by hmmscan"

blastdbcmd \
	-db $answerDB \
	-dbtype prot \
	-entry_batch ${f%%.*}_${jobid}_TMP_Nyeidlist.log \
	-outfmt "%f" \
	| sed 's/lcl|//g' \
	| sed 's/gi|//g' \
	| sed 's/ unnamed protein product//g' \
	> ${f%%.*}_${jobid}_seqs.fasta

numhits=$(grep -c ">" ${f%%.*}_${jobid}_seqs.fasta)
numdoms=$(grep -c ">" ${f%%.*}_${jobid}_domains.fasta)

echo -e "$numhits seqs extracted in general fasta"
echo -e "$numdoms seqs extracted in domains fasta"


echo -e "Taxonomic report..."

while read si ; do
	echo "$si $(grep -c ">${si}_" ${f%%.*}_${jobid}_seqs.fasta)"
done < $sl > ${f%%.*}_${jobid}_seqs.taxonreport

echo -e "\nMultiple domain report...
($(wc -l ${f%%.*}_${jobid}_domains.multiplereport))"
head -n 20 ${f%%.*}_${jobid}_domains.multiplereport
echo "..."

# Pfamsearch (opcional)

if [ $answerPfam = y ] ;
	then echo -e "Computing domain architectures with Pfamscan..."
		perl /home/xavi/Programes/PfamScan/pfam_scan.pl \
		-fasta ${f%%.*}_${jobid}_seqs.fasta \
		-dir /home/xavi/Programes/PfamScan/ \
		> ${f%%.*}_${jobid}_seqs.pfamscan
		while read seqid ; do 
			grep "$seqid" ${f%%.*}_${jobid}_seqs.pfamscan \
			| tr -s  ' ' '\t' \
			| cut -f 7 \
			| xargs
		done < ${f%%.*}_${jobid}_TMP_Nyeidlist.log > ${f%%.*}_${jobid}_TMP_PFAMTABLE_RESULTS
		paste ${f%%.*}_${jobid}_TMP_Nyeidlist.log ${f%%.*}_${jobid}_TMP_PFAMTABLE_RESULTS > ${f%%.*}_${jobid}_seqs.pfamscanarqs
	else echo -e "Skiping Pfamscan search..."
fi

#Ends... 
echo -e "\n\nRemoving temporary files..."
rm *_TMP_*

echo -e "\n\n>> DONE! <<\n\n"

exit

