# input variables
idx="../../data-evolclust/validation_species_map.csv"
faa="../histone_classification/euk.Histone.seqs.fasta"
gen="/home/xavi/Documents/Lab/prova-genomes/"

# input: species
if [ -z $1 ] ; then
	sps_list=$(awk 'NR>1 { print $1 }' ${idx})
else
	sps_list=$1
fi

sps_list="Crei Vcar Caluen Cvar Cocsub Spun Bden Ttra"

mkdir -p data-evolclust/

# loop through species
for sps in $sps_list; do

	timestamp=$(date +%s)

	# define variables for this species
	sph=$(awk '$1 == "'$sps'" { print $2 }' ${idx})
	spa=$(awk '$1 == "'$sps'" { print $1 }' ${idx})

	if [ -s ${gen}/${spa}_long.annot.gtf ] ; then

		# start
		echo "# Characterise candidate genes in ${sph} (syn=${spa})"

		${gen}/${spa}_long.annot.gtf

	else
		echo "# Cannot analyse genes in ${sph} (syn=${spa}); no GTF available"
	fi

done

exit 0

