# inputs
qsub_alignment="qsub_alignment-single.sh" # qsub script to submit alignment & tree jobs
# create dir for searches
ouf="results_archaea/trees_adhoc"
inf="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/"
mkdir -p ${ouf}

# domains list (families with viral hits only)
domains_list="SET DOT1 SIR2 Hist_deacetyl MOZ_SAS Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC SNF2_N"
#domains_list="Acetyltransf_1"

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

# for each of the families with viral hits, get sequence sets (remove sequence redundancy!)
for i in $domains_list ; do

	# fetch virus sequences
	cd-hit -i ${inf}/arc.${i}.domains.fasta -o ${ouf}/arc.${i}.domains.fasta -c 1 -d 0 1> /dev/null
	sed -i "s/>/>arc_/" ${ouf}/arc.${i}.domains.fasta
	sed -i "s/>/>arc_/" ${ouf}/arc.${i}.domains.fasta.clstr
	echo "# cluster $i viral | $(grep -c '>' ${ouf}/arc.${i}.domains.fasta) seqs"

	# blast virus against other domains
	echo "# diamond $i viral to cellular"
	cd-hit -i ${inf}/bac.${i}.domains.fasta -o tmp.bac.fasta -c 0.7 1> /dev/null
	cat ${inf}/euk.${i}.domains.fasta <(sed "s/>/>bac_/" tmp.bac.fasta) > ${ouf}/tmp.${i}.database.fasta
	rm tmp.bac.fasta*
	diamond blastp -q ${ouf}/arc.${i}.domains.fasta -d ${ouf}/tmp.${i}.database.fasta -o ${ouf}/arc.${i}.diamond.to_cellular.csv --more-sensitive --quiet -k 100 --threads 10
	awk '$11 < 1e-5 { print $2 }'  ${ouf}/arc.${i}.diamond.to_cellular.csv | sort -u > ${ouf}/arc.${i}.diamond.to_cellular.txt
        echo "# diamond $i viral to cellular | $(wc -l ${ouf}/arc.${i}.diamond.to_cellular.txt) hits"
	# retrieve sequences from original fasta
	cp ${ouf}/arc.${i}.domains.fasta ${ouf}/arc.${i}.diamond.to_cellular.fasta
	xargs faidx ${ouf}/tmp.${i}.database.fasta < ${ouf}/arc.${i}.diamond.to_cellular.txt  >> ${ouf}/arc.${i}.diamond.to_cellular.fasta

	# do homology clusters
	echo "# mclclus $i"
	do_homology_clusters ${ouf}/arc.${i}.diamond.to_cellular.fasta ${ouf}/arc.${i}.diamond.to_cellular 4 1.15 5
	clu_list=$(awk '{ print $1 }' ${ouf}/arc.${i}.diamond.to_cellular.mcl_filtered.txt | sort -u -V)

	echo "# mclclus $i | launch trees"
	for clu in $clu_list; do
		awk -v clu="$clu" '$1 == clu { print $2 }' ${ouf}/arc.${i}.diamond.to_cellular.mcl_filtered.txt > ${ouf}/arc.${i}.phy.${clu}.seqs.list
		xargs faidx ${ouf}/arc.${i}.diamond.to_cellular.fasta < ${ouf}/arc.${i}.phy.${clu}.seqs.list  > ${ouf}/arc.${i}.phy.${clu}.seqs.fasta
                if [ $(grep -c '>arc_' ${ouf}/arc.${i}.phy.${clu}.seqs.fasta) -gt 0 ] ; then
			qsub -N arc.${i}.phy.${clu} -pe smp 4 $qsub_alignment ${ouf}/arc.${i}.phy.${clu}.seqs.fasta 4 2> /dev/null 1> /dev/null
			echo "# mclclus $i | launch trees | $(grep -c '>' ${ouf}/arc.${i}.phy.${clu}.seqs.fasta) sequences ($(grep -c '>arc_' ${ouf}/arc.${i}.phy.${clu}.seqs.fasta) viral)"
		fi
	done

done
rm ${ouf}/*.clstr

echo "Done"
