# input variables
idx="../data/validation_species.csv" # where to find PFAM-A.hmm database (preformatted)
faa="results_TEfusions/all_hits.fasta"
gli="results_TEfusions/architectures_with_TEs.genes.txt"
pfam="/users/asebe/xgraubove/data/pfam/"
gen="/users/asebe/xgraubove/genomes/data/"
trx="/no_backup/asebe/xgraubove/sra_transcriptomes/"
nth=8

# input: species
# input: species
if [ -z $1 ] ; then
	sps_list="Drer Chabra Rhidel Bralan Chocri Rirr Morvir Spur Emue Notgen Cvel Plafun Hvul Tetwil Sros Aaur Klenit Coka Porpur Sphfal Tcas Crei Exapal Scil Smoe Symmic Acoe Cavfas Hsap Mlei Adig Nvec Bremin Atha"
	# sps_list="Drer Chabra Chocri"
else
	sps_list=$1
fi

mkdir -p results_TEfusions/data

# loop through species
for sps in $sps_list; do

	timestamp=$(date +%s)

	# define variables for this species
	sra_list=$(awk '$3 == "'$sps'" { print $5 }' ${idx} | tr ',' ' ' )
	sph=$(awk '$3 == "'$sps'" { print $1 }' ${idx})
	spa=$(awk '$3 == "'$sps'" { print $3 }' ${idx})
	spa=$(awk '$3 == "'$sps'" { print $3 }' ${idx})
	gen_list=$(grep "^${sph}_"  ${gli} | sort -u)

	# start
	echo "# Characterise candidate fusion genes in ${sph} (syn=${spa})"

	# extract genes of interest
	echo ${gen_list} | xargs faidx -d ' ' ${faa} > results_TEfusions/data/fusions.${sph}.fasta
	echo "# Evaluating $(grep -c ">" results_TEfusions/data/fusions.${sph}.fasta) candidate fusion genes in ${sph}"

	# Pfam scan proper
	if [ -s results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.fasta ] ; then
	echo "# pfamscan already done..."
	else
	echo "# pfamscan..."
	pfam_scan.pl \
		-fasta results_TEfusions/data/fusions.${sph}.fasta \
		-dir ${pfam} \
		-alig \
		-as \
		-cpu 1 > results_TEfusions/data/fusions.${sph}.Pfamscan.csv

	grep -v "^#" results_TEfusions/data/fusions.${sph}.Pfamscan.csv \
		| sed "/^$/d" \
		| tr -s ' ' '\t' \
		| cut -f 1,7 \
		| awk '
		{line="";for (i = 2; i <= NF; i++)
		line = line $i " "; table[$1]=table[$1] line;
		} END {
		for (key in table) print key, table[key];}' \
		| sed "s/ /\t/" \
		| sed "s/ $//" > results_TEfusions/data/fusions.${sph}.Pfamscan.arqdom.csv

	grep -v "^#" results_TEfusions/data/fusions.${sph}.Pfamscan.csv \
		| sed "/^$/d" \
		| tr -s ' ' ' ' \
		| cut -f 1,2,3,6,7 -d ' ' \
		| awk '{
			split($4, a, ".")
			print $1,$2,$3,substr(a[1], 1),$5
			}' \
		> results_TEfusions/data/fusions.${sph}.Pfamscan_TMP_1

	grep "#SEQ" results_TEfusions/data/fusions.${sph}.Pfamscan.csv \
		| sed "/^$/d" \
		| tr -s ' ' '\t' \
		| cut -f 2 \
		| sed "s/-//g" > results_TEfusions/data/fusions.${sph}.Pfamscan_TMP_2

	paste results_TEfusions/data/fusions.${sph}.Pfamscan_TMP_1 results_TEfusions/data/fusions.${sph}.Pfamscan_TMP_2 > results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.csv
	rm -f results_TEfusions/data/fusions.${sph}.*_TMP_*

	awk '{ print ">"$1"|"$2"|"$3"|"$5"\n"$6 }' results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.csv > results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.fasta
	fi

	# find best hit in genome-informed dataset
	echo "# map proteins to ${gen}/${spa}_long.cds.fasta transcripts..."
	diamond blastp -d ${gen}/${spa}_long.pep.fasta --quiet --more-sensitive -q results_TEfusions/data/fusions.${sph}.fasta -o results_TEfusions/data/fusions.${sph}.diamond.csv.tmp
	awk '$3 == 100' results_TEfusions/data/fusions.${sph}.diamond.csv.tmp > results_TEfusions/data/fusions.${sph}.diamond.csv \
	&& rm results_TEfusions/data/fusions.${sph}.diamond.csv.tmp

	echo "# find gene coordinates in ${gen}/${spa}_long.annot.gtf..."
	fgrep -w -f <(awk '$3 == 100 { print $2 }' results_TEfusions/data/fusions.${sph}.diamond.csv) ${gen}/${spa}_long.annot.gtf | awk '$3 == "transcript" { print $1,$4-1,$5-1,$10 }' | tr ' ' '\t' \
	| tr -d '";' | sort -k1,1 -k2,3n > results_TEfusions/data/fusions.${sph}.diamond.annot.bed

	# check gene structure
	echo "# check gene structure in ${gen}/${spa}_long.annot.gtf"
	fgrep -w -f <(awk '$3 == 100 { print $2 }' results_TEfusions/data/fusions.${sph}.diamond.csv) ${gen}/${spa}_long.annot.gtf | bioawk -c gff '$3 == "exon" { print $9 }'  \
	| cut -f1 -d ';'| sed "s/transcript_id \"//" | tr -d '"'| sort | uniq -c | awk '{ print $2,$1 }'| tr ' ' '\t' \
	> results_TEfusions/data/fusions.${sph}.exons_per_gene.csv

	# find NNNN stretches
	# bedtools getfasta...
	# drop genes if there are NNN stretches ANYWHERE
	echo "# gene contiguity in ${gen}/${spa}_gDNA.fasta"
	bedtools getfasta \
		-fi ${gen}/${spa}_gDNA.fasta \
		-bed results_TEfusions/data/fusions.${sph}.diamond.annot.bed \
		-fo  - \
		-name \
	| bioawk -c fastx '{ print $1,$2 }' | grep -v "[nN]+" | awk '{ print $1 }' | sed "s/::.*//" | sort -u \
	| awk 'NR==FNR { l[$2]=$1;next} { for(i = 1; i <= NF; i++) { if ($i in l) $i=l[$i] ; } } { print $0 }' \
		<(awk '$3 == 100 { print $1,$2 }' results_TEfusions/data/fusions.${sph}.diamond.csv) - | sort -u \
	> results_TEfusions/data/fusions.${sph}.nonNstretch.txt

	# find gene clusters in genome assembly (at fixed distance)
	echo "# find gene clusters in ${gen}/${spa}_long.annot.gtf"
	bedtools cluster -i results_TEfusions/data/fusions.${sph}.diamond.annot.bed -d 100000 \
	> results_TEfusions/data/fusions.${sph}.diamond.annot.bed.tmp \
	&& mv results_TEfusions/data/fusions.${sph}.diamond.annot.bed.tmp results_TEfusions/data/fusions.${sph}.diamond.annot.bed

	# check expression of the equivalent transcripts in the salmon quant files
	# more than one transcript will be reported for each protein model
	# record expression status in each sample
	if [ -s results_TEfusions/data/fusions.${sph}.salmon_quant.csv ] ; then
	echo "# check expression already done..."
	else
	echo "# check expression in ${trx}"
	echo -e "Transcript\tLength\tEffectiveLength\tTPM\tNumReads\tSample\tGene" > results_TEfusions/data/fusions.${sph}.salmon_quant.csv
	if [ -s results_TEfusions/data/fusions.${sph}.diamond.csv ] ; then
		for g in ${gen_list} ; do
			gtlist=$(awk '$1 == "'${g}'" { print $2 } ' results_TEfusions/data/fusions.${sph}.diamond.csv)
			for gt in ${gtlist} ; do
			for ex in ${trx}/${spa}/*_salmon/quant.sf ; do
				exc=$(basename $(dirname ${ex}) | cut -f1 -d '_')
				awk '$1 == "'${gt}'" { print $0"\t'${exc}'\t'${g}'" }' ${ex}
			done
			done
		done >> results_TEfusions/data/fusions.${sph}.salmon_quant.csv
	fi
	fi




	# check continuous coverage
	# cross coordinates of domain along CDS sequence (bed file required) with bedgraph coverage to identify regions with drops in coverage
	if [ -s results_TEfusions/data/fusions.${sph}.salmon_cover.csv ] ; then
	echo "# check coverage already done..."
	else
	echo "# check coverage in ${trx}"
	echo -e "Transcript\tstart_cds\tend_cds\tstart_pep\tend_pep\tNumReads" > results_TEfusions/data/fusions.${sph}.salmon_cover.csv
	# find regions with expression
	if [ -s results_TEfusions/data/fusions.${sph}.diamond.csv ] ; then
		cat ${trx}/${spa}/*_salmon.bdg \
		| fgrep -w -f <(awk '$3 == 100 { print $2 }' results_TEfusions/data/fusions.${sph}.diamond.csv) - \
		| sort -k1,1 -k5,5 -k2,2n \
		| awk '$4 > 0' \
		| sort -k1,1 -k2,2n \
		| bedtools merge -i - -c 4 -o mean | awk -v OFS="\t" '{ print $1,$2,$3,int($2/3),int($3/3),$4 }' >> results_TEfusions/data/fusions.${sph}.salmon_cover.csv
	fi
	fi

	# # find regions with zero expression in these transcripts
	# # includes all regions that have zero expression in ALL samples.
	# # excludes regions that have at least one read in at least one sample
	# echo -e "Transcript\tstart_cds\tend_cds\tstart_pep\tend_pep\tNumReads" > results_TEfusions/data/fusions.${sph}.salmon_zerocov.csv
	# bedtools complement \
	# 	-i results_TEfusions/data/fusions.${sph}.salmon_cover.COVER \
	# 	-g <(bioawk -c fastx '{ print $1,length($2) }' ${trx}/${spa}/${spa}_long.cds.fasta | sort -k1,1) \
	# 	-L \
	# 	| awk -v OFS="\t" '{ print $1,$2,$3,int($2/3),int($3/3),0 }' >> results_TEfusions/data/fusions.${sph}.salmon_zerocov.csv && rm results_TEfusions/data/fusions.${sph}.salmon_cover.COVER

	echo "# candidate fusions     = $(grep -c ">" results_TEfusions/data/fusions.${sph}.fasta)"
	echo "# expressed (>0 sample) = $(awk 'NR>1 && $4 > 0 { print $7 }' results_TEfusions/data/fusions.${sph}.salmon_quant.csv | sort -u | wc -l)"
	echo "# contiguous assembly   = $(wc -l results_TEfusions/data/fusions.${sph}.nonNstretch.txt)"
	echo "# ${sps} done!"


	# # blat to genome assembly: map domains to genome assembly (restrict to domains that overlap with best hits in this assembly)
	# if [ -s results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed ] ; then
	# echo "# blat already done..."
	# else
	# echo "# blat to ${gen}/${spa}_gDNA.fasta..."
	# blat ${gen}/${spa}_gDNA.fasta results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.fasta -q=prot -t=dnax -out=blast8 results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.csv
	# awk '$3 >= 99 && $9 < $10 { print $2,$9-1,$10-1,$1 }' results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.csv | tr ' ' '\t' | sort -k1,1 -k2,3n > results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed.tmp
	# bedtools intersect -a results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed.tmp -b results_TEfusions/data/fusions.${sph}.diamond.annot.bed > results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed \
	# && rm results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed.tmp
	# fi
	# echo "# $(cut -f4 results_TEfusions/data/fusions.${sph}.Pfamscan.seqs.blat.bed | cut -f1 -d '|' | sort -u | wc -l) proteins successfully mapped to ${gen}/${spa}_gDNA.fasta..."


done

exit 0

