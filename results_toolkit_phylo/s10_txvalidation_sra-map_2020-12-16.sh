# input variables
idx="../data/validation_species.csv" # where to find PFAM-A.hmm database (preformatted)
out="/no_backup/asebe/xgraubove/sra_transcriptomes"
gen="/users/asebe/xgraubove/genomes/data"
nth=8

# input: species
if [ -z $1 ] ; then
	sps_list="Drer Chocri Morvir Spur Bralan Emue Chabra Rhidel Notgen Cvel Rirr Plafun Hvul Tetwil Sros Aaur Klenit Coka Porpur Sphfal Tcas Crei Exapal Scil Smoe Symmic Acoe Cavfas Hsap Mlei Adig Nvec Bremin"
else
	sps_list=$1
fi

for sps in $sps_list; do

	# define variables for this species
	sra_list=$(awk '$3 == "'$sps'" { print $5 }' ${idx} | tr ',' ' ' )

	# index transcriptome
	mkdir -p ${out}/${sps} && cp ${gen}/${sps}_long.cds.fasta ${out}/${sps}/
	fas=${out}/${sps}/${sps}_long.cds.fasta

	echo "# Index ${fas}"
	# bwa index ${fas}
	if [ -s ${out}/${sps}/salmon_index/versionInfo.json ] ; then
		echo "# Index already present in ${out}/${sps}/salmon_index/"
	else
		echo "# Do index in ${out}/${sps}/salmon_index/"
		salmon index -t ${fas} -i ${out}/${sps}/salmon_index  -k 31 2> /dev/null 1> /dev/null
	fi

	echo "# Launch jobs ${fas}"
	for ri in ${sra_list}; do

		# arguments: run id, sps id, path to fasta (txs), path to output folder, num threads
		qsub \
			-pe smp ${nth} \
			-N "map-${sps}-${ri}" \
			-o ${out}/${sps}/${ri}.log \
			-e ${out}/${sps}/${ri}.log \
			qsub_srasalmon.sh \
				${ri} \
				${sps} \
				${fas} \
				${out}/${sps} \
				${nth}
	done

	echo "# ${sps} done!"

done

exit 0

