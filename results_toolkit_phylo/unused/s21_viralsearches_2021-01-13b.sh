# inputs
qsub_alignment="qsub_alignment-fast.sh" # qsub script to submit alignment & tree jobs
# create dir for searches
ouf="results_viruses/trees_vireuk"
inf="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/"
mkdir -p ${ouf}


# domains list (families with viral hits only)
domains_list="Acetyltransf_1 AF9 BIR Bromodomain Chromo CupinJmjC DOT1 GNAT_acetyltr_2 HIF-1 Hist_deacetyl Histone LinkerHistone PHD PTIP SAM SET SIR2 SNF2_N TUDOR zf-CCHH zf-CXXC"


# for each of the families with viral hits, get sequence sets (remove sequence redundancy!)
for i in $domains_list ; do

	echo "# create dataset $i"
	cat ${inf}/euk.${i}.domains.fasta <(sed "s/>/>vir_/" ${inf}/vir.${i}.domains.fasta) > ${ouf}/cat.${i}.domains.fasta

	# maybe launch whole-dataset quick alignment?
	qsub -N cat.${i} -pe smp 4 -q long-sl7,mem_512_12h -l virtual_free=30G,h_rt=12:00:00 -pe smp 4 qsub_alignment-fast.sh ${ouf}/cat.${i}.domains.fasta 4
	echo "# create dataset $i | $(grep -c ">" ${ouf}/cat.${i}.domains.fasta)"

done

echo "Done"
