# input variables
indx_fn="results_TEfusions/validation_species.csv" # where to find PFAM-A.hmm database (preformatted)
rnad_fn="/home/xavi/dades/data_SRA/"
geno_fn="/home/xavi/dades/genomes/"
nth=8

# mkdirs
reaf="${rnad_fn}"
mkdir -p ${reaf}

function downaload_sra {

	local srr=$1
	local out=$2

	if [ -s "${out}/${srr}_2.fastq.gz" -a -s "${out}/${srr}_2.fastq.gz" ] ; then 
		echo "# Already found: ${out}/${srr}_1.fastq.gz"
		echo "# Already found: ${out}/${srr}_2.fastq.gz"
	else
		# check if file already exists and size is non-zero
		fasterq-dump ${srr} -O ${out} --temp ${out}
		gzip ${out}/${srr}_1.fastq
		gzip ${out}/${srr}_2.fastq 
	fi


}

function map_paired {

	local fas=$1
	local re1=$2
	local re2=$3
	local out=$4
	local pid=$5
	local nth=$6

	# check if index exists, build if absent
	if [ -s "${fas}.sa" ] ; then 
		echo "# Already done: ${fas} index"
	else
		echo "# Do: ${fas} index"
		bwa index ${fas}
	fi

	# check output exists, map if absent
	if [ -s "${out}/${pid}_bwa.f.bam" ] ; then 
		echo "# Already done: ${out}/${pid}_bwa.bam alignments and QC"
	else
		# align
		echo "# Do: ${out}/${pid}_bwa.bam mapping"
		bwa mem \
			-M -t ${nth} \
			${fas} \
			${re1} \
			${re2} \
			-v 0 \
			-R "@RG\tID:${pid}\tSM:${pid}\tPU:nye\tPL:nye\tLB:${pid}" \
			| bamsormadup inputformat=sam threads=${nth} tmpfile=${out}/tmp_$(date +%s) SO=coordinate \
			> ${out}/${pid}_bwa.bam.uf 2> ${out}/${pid}_bwa.err
		samtools index ${out}/${pid}_bwa.bam.uf

		# remove low-quality mapping
		echo "# Do: ${out}/${pid}_bwa.bam QC"
		alignmentSieve \
			--bam ${out}/${pid}_bwa.bam.uf \
			--outFile ${out}/${pid}_bwa.bam.u \
			--numberOfProcessors ${nth} \
			--minMappingQuality 30
		# sort
		bamsormadup inputformat=bam threads=${nth} tmpfile=${out}/tmp_$(date +%s) SO=coordinate < ${out}/${pid}_bwa.bam.u > ${out}/${pid}_bwa.bam \
		&& rm ${out}/${pid}_bwa.bam.u ${out}/${pid}_bwa.bam.uf ${out}/${pid}_bwa.bam.uf.bai

	fi

}

while read -a si ; do
	
	# define variables for this species
	sps=${si[2]}
	fas=${geno_fn}/${sps}_long.cds.fasta
	sra_list=$(echo ${si[4]} | tr ',' ' ' )

	# output folder
	mkdir -p ${reaf}/${sps}

	for ri in ${sra_list}; do

		# first, download reads if needed
		echo "# Download ${sps} $ri"
		downaload_sra ${ri} ${reaf}/${sps}

		# second, map reads to spliced transcripts
		echo "# Map ${sps} ${ri}"
		map_paired \
			${fas} \
			${reaf}/${sps}/${ri}_1.fastq.gz \
			${reaf}/${sps}/${ri}_2.fastq.gz \
			${reaf}/${sps}/ \
			${ri} \
			${nth}

		# third, extract coverage along each transcript
		echo "# Coverage ${sps} ${ri}"
		bedtools genomecov -ibam ${reaf}/${sps}/${ri}_bwa.bam > ${reaf}/${sps}/${ri}_bwa.bdg

	done
	
done < <( awk 'NR > 1 && $5 != "NA" && ($3 == "Drer" || $3 == "Emue" || $3 == "Chocri" || $3 == "Morvir")' ${indx_fn} )

echo "# All done!"
exit 0