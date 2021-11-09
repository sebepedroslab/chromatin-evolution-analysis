#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7,mem_512_12h
#$ -l virtual_free=20G,h_rt=43200
#$ -o tmp/
#$ -e tmp/

ri=$1   # run id (from SRA)
sps=$2  # sps id
fas=$3  # fasta to map
out=$4  # output folder path
nth=$5  # num threads

# function
function downaload_sra {

	local srr=$1
	local out=$2
	local nth=$3

	if [ -s "${out}/${srr}_2.fastq.gz" -a -s "${out}/${srr}_2.fastq.gz" ] ; then 
		echo "# Already found: ${out}/${srr}_1.fastq.gz"
		echo "# Already found: ${out}/${srr}_2.fastq.gz"
	else
		# check if file already exists and size is non-zero
		fasterq-dump ${srr} -O ${out} --temp ${out} -e ${nth}
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
			| /users/asebe/xgraubove//miniconda3/envs/atacpip/bin/bamsort \
			inputformat=sam threads=${nth} tmpfile=${out}/tmp_$(date +%s) SO=coordinate \
			> ${out}/${pid}_bwa.bam.uf 2> ${out}/${pid}_bwa.err
		samtools index ${out}/${pid}_bwa.bam.uf

		# remove low-quality mapping
		echo "# Do: ${out}/${pid}_bwa.bam QC"
		/users/asebe/xgraubove//miniconda3/envs/atacpip/bin/alignmentSieve \
			--bam ${out}/${pid}_bwa.bam.uf \
			--outFile ${out}/${pid}_bwa.bam.u \
			--numberOfProcessors ${nth} \
			--minMappingQuality 30
		# sort
		/users/asebe/xgraubove//miniconda3/envs/atacpip/bin/bamsort \
		inputformat=bam threads=${nth} tmpfile=${out}/tmp_$(date +%s) SO=coordinate < ${out}/${pid}_bwa.bam.u > ${out}/${pid}_bwa.bam \
		&& rm ${out}/${pid}_bwa.bam.u ${out}/${pid}_bwa.bam.uf ${out}/${pid}_bwa.bam.uf.bai

	fi

}

# first, download reads if needed
echo "# Download ${sps} $ri"
mkdir -p ${out}
downaload_sra ${ri} ${out} ${nth}

# second, map reads to spliced transcripts
echo "# Map ${sps} ${ri}"
map_paired \
	${fas} \
	${out}/${ri}_1.fastq.gz \
	${out}/${ri}_2.fastq.gz \
	${out}/ \
	${ri} \
	${nth}

echo "# Coverage ${sps} ${ri}"
bedtools genomecov -ibam ${out}/${ri}_bwa.bam > ${out}/${ri}_bwa.bdg