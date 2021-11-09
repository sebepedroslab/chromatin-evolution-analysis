#!/bin/bash

###Requirements: Install https://pypi.python.org/pypi/bioscripts.convert/0.4
#				 Install Pfamscan.pl
#				 Add bin/hmmer-3.2.1/easel/miniapps to your path.


# Runs in a folder with .hmm profiles
# First stdin is 00_taxon_index with taxon abbreviations
# Second stdin is domain extraction mode: "all" or "single" (optional; default is all)
# Third stdin is number of flanking aminoacid to fetch with the domain(optional; default is 0)
# Fourth stdin is E-value in hmmer suitable format (optional; default is 0.00001)


rm -f *_TMP_*

#Defining arguments
if [ -z "$1" ]
	then echo "\nERROR! You need to specify arguments!First stdin is 00_taxon_index with taxon abbreviations\nSecond stdin is E-value in hmmer suitable format (optional; default is 0.00001)\n" ; exit
fi

echo "The taxonomic index database is $1" ; TaxonIndex="$1"

if [ -z "$2" ]
	then Domain_mode="$2"
	else Domain_mode="all"
fi


if ! [[ "$Domain_mode" =~ ("single"|"all") ]]
	then echo "Valid domain extraction modes are 'all' or 'single'"
	exit
else 
	echo "I will extract $Domain_mode domains"
fi

if [ -z "$3" ]
	then Flank=0
	else Flank="$3"
fi

if [ -z "$4" ]
	then echo "The default hmmsearch uses Pfam gathering threshold (--cut_ga)" ;
	else echo "The custom hmmsearch E-value is $4" ; Eval="$4"
fi


echo "Extracted domain hits will be extended by $Flank aminoacids"

echo "\nStandard databases to search:
	Eukaryotes -> /users/asebe/asebe/proj/other/euk_proteomes_last/blast_db/all_euks2.fasta
	Archaea -> /Users/arnau/Documents/DataBases/Proteomes_Archaea/AllArcahea.fasta 
	Bacteria -> /Users/arnau/Documents/DataBases/Proteomes_Bacteria/BacteriaNonDraft.fasta "


echo "\nThen, give me your protein database (fasta with BlastDB format):"
read answerDB

cut -f 2 $TaxonIndex > _TMP_TaxonIndex2

#Pfamscan for hits (computes prot architecture)
echo "\nDo you want to run domain architecture analyses with Pfamscan? It takes a while (y/n)"
read answerPfam
echo "\n\n"

#Look for hmmsearch hits, extract hits, computes architectures, computes taxonomic report

	for f in *.hmm ; do 
	
		echo "\n>>> Looking for $f with hmmsearch in $answerDB <<<\n"
		
		if [ -z $Eval ]; then 
			hmmsearch --domtblout ${f%%.*}_out.dtbl  -A ${f%%.*}_TMP_align.stockholm --cpu 16 --cut_ga $f $answerDB > ${f%%.*}.hmmerout
		else
			hmmsearch --domtblout ${f%%.*}_out.dtbl  -A ${f%%.*}_TMP_align.stockholm --cpu 16 --domE $Eval $f $answerDB > ${f%%.*}.hmmerout
		fi
		
		grep subseq ${f%%.*}_TMP_align.stockholm|tr -s " " "\t"|cut -f 6|sort -u > ${f%%.*}_TMP_Nyeidlist.log

		numseqs=`sort -u ${f%%.*}_TMP_Nyeidlist.log | sed '/^$/d' | wc -l | sed -e 's/^[ \t]*//'`


		echo "Extracting hits..."
		blastdbcmd -db $answerDB -dbtype prot -entry_batch ${f%%.*}_TMP_Nyeidlist.log -outfmt "%f" | sed 's/lcl|//g' | sed 's/ unnamed protein product//g' |tr -d ' ' > ${f%%.*}_${numseqs}seqs.fasta
		echo "Taxonomic report..."
		
		numhits=`grep -c ">" ${f%%.*}_${numseqs}seqs.fasta`
		if [ $numseqs == $numhits ] 
			then echo "$numseqs hits found..."
			else echo "\n\n>> ERROR, I found $numseqs seqs and $numhits hits and cannot extract all of them! <<\n\n" ; exit
		fi
		
		echo "Fetching domains..."
		
		if [ "$Domain_mode" == "single" ]; then
			sort -r -g -k1,16 ${f%%.*}_out.dtbl > ${f%%.*}_out2.dtbl ##VERY IMPORTANT: we neet to reorder the hmmeroutput to pick the domain with the lowest Evalue when >1 domain per prot (if single mode)
			grep -v "^#" ${f%%.*}_out2.dtbl | awk -v flank="$Flank" '{print $1, $20-flank, $21+flank, $1}'|awk '$2<1 {$2=1}1' > ${f%%.*}_TMP_coordinates
			cut -f 1 ${f%%.*}_TMP_coordinates -d ' '|sort|uniq -d > ${f%%.*}_CAREFUL_prots_with_multiple_domains
			while read id; do grep -w -m 1 $id ${f%%.*}_TMP_coordinates >> ${f%%.*}_TMP_coordinates2;done < ${f%%.*}_TMP_Nyeidlist.log
		elif [ "$Domain_mode" == "all" ]; then
			sort -g -k1,12 ${f%%.*}_out.dtbl > ${f%%.*}_out2.dtbl   #in this case we sort by coordinates
			grep -v "^#" ${f%%.*}_out2.dtbl | awk -v flank="$Flank" '{print $1, $20-flank, $21+flank, $1}'|awk '$2<1 {$2=1}1' > ${f%%.*}_TMP_coordinates
			cut -f 1 ${f%%.*}_TMP_coordinates -d ' '|sort|uniq -d > ${f%%.*}_CAREFUL_prots_with_multiple_domains
			while read id; do grep -w $id ${f%%.*}_TMP_coordinates | awk '{print $1"\t"$2"\t"$3}' >> ${f%%.*}_TMP_coordinates2_mod_name;done < ${f%%.*}_TMP_Nyeidlist.log #awk because we need to modify only the 1st column ID
			
			while read dup_id; do perl -pi -e 's/\b'${dup_id}'\b/$&.'_domain'.++$A /ge' ${f%%.*}_TMP_coordinates2_mod_name;done < ${f%%.*}_CAREFUL_prots_with_multiple_domains  #append a domain_number label to ids with > 1 domain
			
			paste ${f%%.*}_TMP_coordinates2_mod_name ${f%%.*}_TMP_coordinates | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${f%%.*}_TMP_coordinates2   ##for esl-sfetch, first column must be the specific domain name and the 4th one the general protein ID
			
		fi
		
		#awk '{print $1"_"$2"_"$3, $2, $3, $1"_"$2"_"$3}' ${f%%.*}_TMP_coordinates2 > ${f%%.*}_TMP_coordinates3; if we wanted to add coordinates to prot name, annoying.
		
		
		bioawk -c fastx '{print $name,length($seq)}' ${f%%.*}_${numseqs}seqs.fasta > ${f%%.*}_TMP_seq_lengths
		while read line;do  prot_id=$(echo $line|cut -f 1 -d " "); grep -w $prot_id ${f%%.*}_TMP_seq_lengths >> ${f%%.*}_TMP_seq_lengths2 ;done < ${f%%.*}_TMP_coordinates
		paste ${f%%.*}_TMP_coordinates2 ${f%%.*}_TMP_seq_lengths2 > ${f%%.*}_TMP_coordinates3
		awk '$3>$6 {$3=$6}1' ${f%%.*}_TMP_coordinates3|awk '{print $1,$2,$3,$4}' > ${f%%.*}_coordinates
		
		esl-sfetch --index ${f%%.*}_${numseqs}seqs.fasta
		esl-sfetch -C -f ${f%%.*}_${numseqs}seqs.fasta ${f%%.*}_coordinates > ${f%%.*}_domains.fasta
		
		
		
		while read taxonab ; do
			grep -c ">$taxonab\_" ${f%%.*}_${numseqs}seqs.fasta 
		done < _TMP_TaxonIndex2 >  ${f%%.*}_TMP_TaxonReportCount.log
		paste _TMP_TaxonIndex2 ${f%%.*}_TMP_TaxonReportCount.log > ${f%%.*}_TMP_TaxonReport.log

		if [ $answerPfam == y ] ;
			then echo "Computing domain architectures with Pfamscan..."
				pfam_scan.pl -fasta ${f%%.*}_${numseqs}seqs.fasta -dir /nfs/users2/asebe/asebe/proj/other/Pfam > ${f%%.*}_Pfamscan.csv
				while read seqid ; do 
					grep "$seqid" ${f%%.*}_Pfamscan.csv | tr -s  ' ' '\t' | cut -f 7 | xargs 
				done < ${f%%.*}_TMP_Nyeidlist.log > ${f%%.*}_TMP_PFAMTABLE_RESULTS
				paste ${f%%.*}_TMP_Nyeidlist.log ${f%%.*}_TMP_PFAMTABLE_RESULTS > ${f%%.*}_DomainArchitectures.csv 
			else echo "Skiping Pfamscan search..."
		fi
	done


#Pastes together taxonomic reports 
echo "\n\nBuilding global taxon report..."

for f in *_TMP_TaxonReport.log ; do 
	echo "${f%%_TMP_TaxonReport.*}" >> ${f%%.*}_TMP_lol 
	cut -f 2 $f >> ${f%%.*}_TMP_lol 
done 

paste *_TMP_lol > _TMP_TaxonIndexNumbers
echo "Genomes" > _TMP_TaxonIndexTitle
cat _TMP_TaxonIndexTitle _TMP_TaxonIndex2 > _TMP_TaxonIndex3 
#if [ $answerStandard == y ] ; then echo "ARCHAEA\nBACTERIA" >> _TMP_TaxonIndex3 ; fi
sed '1d' _TMP_TaxonIndex3 > _TMP_TaxonIndex4
paste _TMP_TaxonIndex3 _TMP_TaxonIndexNumbers | tr '\t' ';' > 00_TaxonReport.csv


#Ends... 
echo "\n\nRemoving temporary files..."
rm *_TMP_*
echo "Removing empty fasta files..."
rm -f *_0seqs_*.fasta


echo "\n\n>> DONE! <<\n\n"

exit

