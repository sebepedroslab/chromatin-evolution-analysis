# inputs
qsub_alignment="qsub_alignment-single.sh" # qsub script to submit alignment & tree jobs
# create dir for searches
ouf="results_preeuk/trees_adhfam"
inf="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/"
mkdir -p ${ouf}

# domains list (families with viral hits only)
#domains_list="SET DOT1 SIR2 Hist_deacetyl MOZ_SAS Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC SNF2_N"
#domains_list="Acetyltransf_1"
domains_list="ASF1_hist_chap DOT1 Hist_deacetyl MOZ_SAS SIR2 SET Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC ING MBT PWWP SNF2_N Histone YEATS Nucleoplasmin"
#domains_list="Nucleoplasmin"

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
	cd-hit -i ${inf}/euk.${i}.domains.fasta -o ${ouf}/euk.${i}.domains.fasta -c 1 -d 0 1> /dev/null
	echo "# cluster $i euk | $(grep -c '>' ${ouf}/euk.${i}.domains.fasta) seqs"

        # blast euks against bacteria and archaea
        echo "# reduce prok $i | cd-hit"
        cd-hit -i ${inf}/bac.${i}.domains.fasta -o tmp.bac.fasta -c 0.9 1> /dev/null
        cd-hit -i ${inf}/arc.${i}.domains.fasta -o tmp.arc.fasta -c 0.95 1> /dev/null
        cat <(sed "s/>/>bac_/" tmp.bac.fasta) <(sed "s/>/>arc_/" tmp.arc.fasta) > ${ouf}/tmp.${i}.database.fasta
        #rm tmp.bac.fasta* tmp.arc.fasta*

	diamond blastp -q ${ouf}/euk.${i}.domains.fasta -d ${ouf}/tmp.${i}.database.fasta -o ${ouf}/euk.${i}.diamond.to_cellular.csv --more-sensitive --quiet -k 10 --threads 10
	awk '$11 < 1e-5 { print $2 }'  ${ouf}/euk.${i}.diamond.to_cellular.csv | sort -u > ${ouf}/euk.${i}.diamond.to_cellular.txt
        echo "# diamond $i euk to prok | $(wc -l ${ouf}/euk.${i}.diamond.to_cellular.txt) hits"

	# retrieve sequences from original fasta
	cp ${ouf}/euk.${i}.domains.fasta ${ouf}/euk.${i}.diamond.to_cellular.fasta
	xargs faidx ${ouf}/tmp.${i}.database.fasta < ${ouf}/euk.${i}.diamond.to_cellular.txt  >> ${ouf}/euk.${i}.diamond.to_cellular.fasta

	# do homology clusters
	echo "# mclclus $i"
	do_homology_clusters ${ouf}/euk.${i}.diamond.to_cellular.fasta ${ouf}/euk.${i}.diamond.to_cellular 4 1.2 5
	clu_list=$(awk '{ print $1 }' ${ouf}/euk.${i}.diamond.to_cellular.mcl_filtered.txt | sort -u -V)

	echo "# mclclus $i | launch trees"
	for clu in $clu_list; do
		awk -v clu="$clu" '$1 == clu { print $2 }' ${ouf}/euk.${i}.diamond.to_cellular.mcl_filtered.txt > ${ouf}/euk.${i}.phy.${clu}.seqs.list
		xargs faidx ${ouf}/euk.${i}.diamond.to_cellular.fasta < ${ouf}/euk.${i}.phy.${clu}.seqs.list  > ${ouf}/euk.${i}.phy.${clu}.seqs.fasta
                if [ $(grep -c '>bac_' ${ouf}/euk.${i}.phy.${clu}.seqs.fasta) -gt 0 ] || [ $(grep -c '>arc_' ${ouf}/euk.${i}.phy.${clu}.seqs.fasta) -gt 0 ] ; then
			qsub -N euk.${i}.phy.${clu} -pe smp 4 $qsub_alignment ${ouf}/euk.${i}.phy.${clu}.seqs.fasta 4 2> /dev/null 1> /dev/null
			echo "# mclclus $i | launch trees | $(grep -c '>' ${ouf}/euk.${i}.phy.${clu}.seqs.fasta) sequences ($(grep -c '>arc_' ${ouf}/euk.${i}.phy.${clu}.seqs.fasta) arc, $(grep -c '>bac_' ${ouf}/euk.${i}.phy.${clu}.seqs.fasta) bac)"
		fi
	done

done
rm ${ouf}/*.clstr

echo "Done"
