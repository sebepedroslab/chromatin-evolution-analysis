# input variables
idx="../../data/validation_species_map.csv"
faa="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/euk.Histone.seqs.fasta"
gen="/users/asebe/xgraubove/genomes/data/"

# input: species
if [ -z $1 ] ; then
	sps_list=$(awk 'NR>1 { print $1 }' ${idx})
else
	sps_list=$1
fi

mkdir -p data

# loop through species
for sps in $sps_list; do

	timestamp=$(date +%s)

	# define variables for this species
	sph=$(awk '$1 == "'$sps'" { print $2 }' ${idx})
	spa=$(awk '$1 == "'$sps'" { print $1 }' ${idx})
	gen_list=$(bioawk -c fastx '{ print $1 }' ${faa} | grep "^${sph}_" | sort -u)

	if [ -s ${gen}/${spa}_long.annot.gtf ] ; then
		# start
		echo "# Characterise candidate genes in ${sph} (syn=${spa})"

		# extract genes of interest
		echo ${gen_list} | xargs faidx -d ' ' ${faa} > data/genes.${sph}.fasta
		echo "# Evaluating $(grep -c ">" data/genes.${sph}.fasta) candidate genes in ${sph}"

		# find best hit in genome-informed dataset
		echo "# map proteins to ${gen}/${spa}_long.cds.fasta transcripts..."
		diamond blastp -d ${gen}/${spa}_long.pep.fasta --quiet --more-sensitive -q data/genes.${sph}.fasta -o data/genes.${sph}.diamond.csv.tmp
		awk '$3 == 100' data/genes.${sph}.diamond.csv.tmp > data/genes.${sph}.diamond.csv \
		&& rm data/genes.${sph}.diamond.csv.tmp

		echo "# find gene coordinates in ${gen}/${spa}_long.annot.gtf..."
		fgrep -w -f <(awk '$3 == 100 { print $2 }' data/genes.${sph}.diamond.csv) ${gen}/${spa}_long.annot.gtf | awk '$3 == "transcript" { print $1,$4-1,$5-1,$10,$7 }' | tr ' ' '\t' \
		| tr -d '";' | sort -k1,1 -k2,2n  | bedtools sort > data/genes.${sph}.diamond.annot.bed

		# find gene clusters in genome assembly (at fixed distance)
		echo "# find gene clusters in ${gen}/${spa}_long.annot.gtf"
		bedtools cluster -i data/genes.${sph}.diamond.annot.bed -d 100000 \
		> data/genes.${sph}.diamond.annot.bed.tmp \
		&& mv data/genes.${sph}.diamond.annot.bed.tmp data/genes.${sph}.diamond.annot.bed

		# find closest genes in genome assembly
		echo "# find closest genes clusters in ${gen}/${spa}_long.annot.gtf"
		bedtools closest -a data/genes.${sph}.diamond.annot.bed -b <(awk '$3 == "transcript"  { print $1,$4-1,$5-1,$10,$7 }' ${gen}/${spa}_long.annot.gtf | tr ' ' '\t' | tr -d '";' | sort -k1,1 -k2,2n | bedtools sort) \
		-io -fd -D ref > data/genes.${sph}.diamond.annot.closest.downstream.bed
		bedtools closest -a data/genes.${sph}.diamond.annot.bed -b <(awk '$3 == "transcript"  { print $1,$4-1,$5-1,$10,$7 }' ${gen}/${spa}_long.annot.gtf | tr ' ' '\t' | tr -d '";' | sort -k1,1 -k2,2n | bedtools sort) \
		-io -fu -D ref > data/genes.${sph}.diamond.annot.closest.upstream.bed

	else
		echo "# Cannot analyse genes in ${sph} (syn=${spa}); no GTF available"
	fi

done

exit 0

