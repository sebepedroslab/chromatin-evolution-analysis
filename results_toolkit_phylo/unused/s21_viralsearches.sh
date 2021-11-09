# inputs
qsub_alignment="qsub_alignment.sh" # qsub script to submit alignment & tree jobs
# create dir for searches
mkdir -p results_viruses/sequences

# domains list (families with viral hits only)
domains_list="Acetyltransf_1 AF9 BIR Bromodomain Chromo CupinJmjC DOT1 GNAT_acetyltr_2 HIF-1 Hist_deacetyl Histone Kelch LinkerHistone PHD PTIP SAM SET SIR2 SNF2_N TUDOR WD40 zf-C2H2 zf-CCHH zf-CXXC"
domains_list="SET BIR Bromodomain Chromo CupinJmjC DOT1 GNAT_acetyltr_2 HIF-1 Histone LinkerHistone PHD PTIP SIR2 SNF2_N TUDOR zf-CCHH zf-CXXC Hist_deacetyl Acetyltransf_1"


# for each of the families with viral hits, get sequence sets (remove sequence redundancy!)
for i in $domains_list ; do

	echo "# create dataset $i"
	cd-hit -i /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/euk.${i}.domains.fasta -o results_viruses/sequences/cat.${i}.domains.fasta -c 0.95 1> /dev/null
	cat /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/vir.${i}.domains.fasta > results_viruses/sequences/vir.${i}.domains.fasta
	sed "s/>/>vir_/" results_viruses/sequences/vir.${i}.domains.fasta >> results_viruses/sequences/cat.${i}.domains.fasta && rm results_viruses/sequences/vir.${i}.domains.fasta
	cd-hit -i /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/arc.${i}.domains.fasta -o results_viruses/sequences/arc.${i}.domains.fasta -c 0.95 1> /dev/null
	sed "s/>/>arc_/" results_viruses/sequences/arc.${i}.domains.fasta >> results_viruses/sequences/cat.${i}.domains.fasta && rm results_viruses/sequences/arc.${i}.domains.fasta
	cd-hit -i /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/bac.${i}.domains.fasta -o results_viruses/sequences/bac.${i}.domains.fasta -c 0.90 1> /dev/null
	sed "s/>/>bac_/" results_viruses/sequences/bac.${i}.domains.fasta >> results_viruses/sequences/cat.${i}.domains.fasta && rm results_viruses/sequences/bac.${i}.domains.fasta

	# maybe launch whole-dataset quick alignment?
	# qsub -N cat.${i} -pe smp 4 -q mem_512_12h,long-sl7 -l virtual_free=30G,h_rt=12:00:00 -pe smp 4 qsub_alignment-fast.sh results_viruses/sequences/cat.${i}.domains.fasta 4
done
rm results_viruses/sequences/*.clstr



# for each family, partition groups
function do_homology_clusters {

    local fas=$1
    local out=$2
    local cpu=$3
    local inf=$4
    local min=$5

    # pairwise alignments
    diamond makedb --in ${fas} -d ${fas} --quiet
    diamond blastp \
        --more-sensitive \
		--max-target-seqs 100 \
        -d ${fas} \
        -q ${fas} \
        -o ${out}.diamond.csv \
        --quiet \
        --threads ${cpu}

    # partition sequences with MCL (edge weights are alignment bitscores)
    # we select a very low inflation value so as to obtain large, non-granular clusters
    awk '{ print $1,$2,$12 }' ${out}.diamond.csv > ${out}.diamond.abc
    mcl ${out}.diamond.abc --abc -I ${inf} -o ${out}.mcl.csv 2> /dev/null

    # compress stuff
    gzip -f ${out}.diamond.abc
    gzip -f ${out}.diamond.csv

    # ignore homology clusters with <X seqs
    awk 'NF > '$min'' ${out}.mcl.csv > ${out}.mcl_filtered.csv

    # remove clusters with <X seqs
    awk '{ for(i = 1; i <= NF; i++) { print "HG"NR,$i } }' ${out}.mcl_filtered.csv > ${out}.mcl_filtered.txt

}

# partition clusters
for i in $domains_list ; do

	echo "# launch search $i | homology groups with viral sequences..."
	do_homology_clusters results_viruses/sequences/cat.${i}.domains.fasta results_viruses/sequences/cat.${i} 4 1.5 5
	clu_list=$(grep " vir_" results_viruses/sequences/cat.${i}.mcl_filtered.txt | awk '{ print $1 }' | sort -u -V)

	for clu in $clu_list; do
		awk -v clu="$clu" '$1 == clu { print $2 }' results_viruses/sequences/cat.${i}.mcl_filtered.txt > results_viruses/sequences/cat.${i}.${clu}.seqs.list
		echo "# launch search $i | $clu $(wc -l results_viruses/sequences/cat.${i}.${clu}.seqs.list)"
		xargs faidx results_viruses/sequences/cat.${i}.domains.fasta  < results_viruses/sequences/cat.${i}.${clu}.seqs.list > results_viruses/sequences/cat.${i}.${clu}.seqs.fasta
		qsub -N vir.${i}.${clu} -q mem_512_12h,long-sl7 -l virtual_free=30G,h_rt=12:00:00 -pe smp 4 $qsub_alignment results_viruses/sequences/cat.${i}.${clu}.seqs.fasta 4
	done

done

echo "Done"
