#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7,mem_512_12h
#$ -l virtual_free=20G,h_rt=43200

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

	if [ -s "${out}/${srr}_2.fastq.gz" ] || [ -s "${out}/${srr}.fastq.gz" ] ; then 
		echo "# Already found: ${out}/${srr} reads $(ls ${out}/${srr}*.fastq.gz)"
	else
		# check if file already exists and size is non-zero
		fasterq-dump ${srr} -O ${out} --temp ${out} -e ${nth}
		if [ -s "${out}/${srr}_1.fastq" ] ; then
			gzip ${out}/${srr}_1.fastq
			gzip ${out}/${srr}_2.fastq 
		else
			gzip ${out}/${srr}.fastq 
		fi
	fi


}

function salmon_paired {

	local fas=$1
	local re1=$2
	local re2=$3
	local out=$4
	local pid=$5
	local nth=$6

	# check if index exists, build if absent
	if [ -s "${out}/salmon_index/versionInfo.json" ] ; then 
		echo "# Already done: ${fas} index"
	else
		echo "# Do: ${fas} index"
		salmon index -t ${fas} -i salmon_index -k 31
	fi

	# check output exists, map if absent
	if [ -s "${out}/${pid}_salmon/quant.sf" ] ; then 
		echo "# Already done: ${out}/${pid}_salmon/quant.sf"
	else
		# align
		echo "# Do: ${out}/${pid}_bwa.bam mapping"
		rm -rf ${out}/${pid}_salmon
		salmon quant -i ${out}/salmon_index -l A -1 ${re1} -2 ${re2}  -o ${out}/${pid}_salmon --writeMappings=${out}/${pid}_salmon.sam -p ${nth}
		/users/asebe/xgraubove/miniconda3/envs/atacpip/bin/bamsort inputformat=sam threads=${nth} tmpfile=${out}/tmp_${pid}_$(date +%s) SO=coordinate < ${out}/${pid}_salmon.sam > ${out}/${pid}_salmon.bam && rm ${out}/${pid}_salmon.sam
	fi

}

function salmon_unpaired {

	local fas=$1
	local re1=$2
	local out=$3
	local pid=$4
	local nth=$5

	# check if index exists, build if absent
	if [ -s "${out}/salmon_index" ] ; then 
		echo "# Already done: ${fas} index"
	else
		echo "# Do: ${fas} index"
		salmon index -t ${fas} -i salmon_index  -k 31
	fi

	# check output exists, map if absent
	if [ -s "${out}/${pid}_salmon/quant.sf" ] && [ -s "${out}/${pid}_salmon.bam" ] ; then 
		echo "# Already done: ${out}/${pid}_salmon/quant.sf"
	else
		# align
		echo "# Do: ${out}/${pid}_bwa.bam mapping"
		rm -rf ${out}/${pid}_salmon
		salmon quant -i ${out}/salmon_index -l A -r ${re1} -o ${out}/${pid}_salmon --writeMappings=${out}/${pid}_salmon.sam -p ${nth}
		/users/asebe/xgraubove/miniconda3/envs/atacpip/bin/bamsort inputformat=sam threads=${nth} tmpfile=${out}/tmp_${pid}_$(date +%s) SO=coordinate < ${out}/${pid}_salmon.sam > ${out}/${pid}_salmon.bam && rm ${out}/${pid}_salmon.sam
	fi

}



# first, download reads if needed
echo "# Download ${sps} $ri"
mkdir -p ${out}
downaload_sra ${ri} ${out} ${nth}

# second, map reads to spliced transcripts
if [ -s "${out}/${ri}_1.fastq.gz" -a -s "${out}/${ri}_2.fastq.gz" ] ; then
	echo "# Map ${sps} ${ri} (paired)"
	salmon_paired   ${fas} ${out}/${ri}_1.fastq.gz ${out}/${ri}_2.fastq.gz ${out}/ ${ri} ${nth}
else
	echo "# Map ${sps} ${ri} (unpaired)"
	salmon_unpaired ${fas} ${out}/${ri}.fastq.gz                           ${out}/ ${ri} ${nth}
fi

echo "# Coverage ${sps} ${ri}"
if [ -s "${out}/${ri}_salmon.bdg" ] ; then
	echo "# Bedgraph found in ${out}/${ri}_salmon.bdg, skip"
else
	bedtools genomecov -bga -ibam ${out}/${ri}_salmon.bam > ${out}/${ri}_salmon.bdg
fi
