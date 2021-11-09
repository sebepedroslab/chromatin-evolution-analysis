# inputs
qsub_alignment="qsub_alignment-single.sh" # qsub script to submit alignment & tree jobs
# create dir for searches
ouf="results_preeuk/trees_adhog"
inf="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/"
mkdir -p ${ouf}

# domains list (families with viral hits only)
#domains_list="SET DOT1 SIR2 Hist_deacetyl MOZ_SAS Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC SNF2_N Histone"
domains_list="ASF1_hist_chap DOT1 Hist_deacetyl MOZ_SAS SIR2 SET Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC SNF2_N Histone ING MBT PWWP"
domains_list="Acetyltransf_1 GNAT_acetyltr_2 CupinJmjC SNF2_N Histone ING MBT PWWP"

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
    awk '{ for(i = 1; i <= NF; i++) { print "SHG"NR,$i } }' ${out}.mcl_filtered.csv > ${out}.mcl_filtered.txt

}

# for each of the families with viral hits, get sequence sets (remove sequence redundancy!)
for i in $domains_list ; do

	# blast euks against bacteria and archaea
        echo "# reduce prok $i | cd-hit"
	cd-hit -i ${inf}/bac.${i}.domains.fasta -o tmp.bac.fasta -c 0.9 1> /dev/null
	cd-hit -i ${inf}/arc.${i}.domains.fasta -o tmp.arc.fasta -c 0.95 1> /dev/null
	cat <(sed "s/>/>bac_/" tmp.bac.fasta) <(sed "s/>/>arc_/" tmp.arc.fasta) > ${ouf}/tmp.${i}.database.fasta
	rm tmp.bac.fasta* tmp.arc.fasta*

	# fetch euk sequences, separated by OG
	list_ogs=$(awk 'NR>1 { print $2 }' orthogroups_euk.csv | grep "${i}.HG" | sort | uniq -c | awk '$1 > 20 { print $2 }'| sort -u)
	for h in ${list_ogs} ; do

		hc=$(echo ${h} | cut -f1-3 -d '.' | cut -f1 -d ':')
		# get eukaryotic proteins in this OG
		grep -w ${h} orthogroups_euk.csv | cut -f1 | sort -u > ${ouf}/euk.${hc}.seed.txt
		xargs faidx ${inf}/euk.${i}.seqs.fasta < ${ouf}/euk.${hc}.seed.txt > ${ouf}/euk.${hc}.seed.fasta
                echo "# get seed $i | ${hc} | $(grep -c '>' ${ouf}/euk.${hc}.seed.fasta ) seed"

		# find similar sequences in other OGs
		bioawk -c fastx '{ print $1, $2 }' ${inf}/euk.${i}.domains.fasta | fgrep -v -f ${ouf}/euk.${hc}.seed.txt | awk '{ print ">"$1"\n"$2 }' > ${ouf}/tmp.euk.fasta
		diamond blastp -q ${ouf}/euk.${hc}.seed.fasta -d ${ouf}/tmp.euk.fasta -o ${ouf}/euk.${hc}.seed.diamond.othereuk.csv --more-sensitive --quiet -k 20 --threads 10
		awk '$11 < 1e-5 { print $2 }' ${ouf}/euk.${hc}.seed.diamond.othereuk.csv | fgrep -v -f ${ouf}/euk.${hc}.seed.txt | sort -u > ${ouf}/euk.${hc}.seed.diamond.othereuk.txt
                xargs faidx ${inf}/euk.${i}.domains.fasta < ${ouf}/euk.${hc}.seed.diamond.othereuk.txt | sed "s/>/>eukoth_/"> ${ouf}/euk.${hc}.seed.diamond.othereuk.fasta
                echo "# get other eukaryotes $i | ${hc} | $(grep -c '>' ${ouf}/euk.${hc}.seed.diamond.othereuk.fasta) hits"

                # find similar sequences in prokaryotes
                diamond blastp -q ${ouf}/euk.${hc}.seed.fasta -d ${ouf}/tmp.${i}.database.fasta -o ${ouf}/euk.${hc}.seed.diamond.prokaryotes.csv --more-sensitive --quiet -k 20 --threads 10
                awk '$11 < 1e-5 { print $2 }' ${ouf}/euk.${hc}.seed.diamond.prokaryotes.csv | sort -u > ${ouf}/euk.${hc}.seed.diamond.prokaryotes.txt
                xargs faidx ${ouf}/tmp.${i}.database.fasta < ${ouf}/euk.${hc}.seed.diamond.prokaryotes.txt > ${ouf}/euk.${hc}.seed.diamond.prokaryotes.fasta
                echo "# get prokaryotes $i | ${hc} | $(grep -c '>' ${ouf}/euk.${hc}.seed.diamond.prokaryotes.fasta) hits"

		# concatenate all hits
		cat  ${ouf}/euk.${hc}.seed.fasta  ${ouf}/euk.${hc}.seed.diamond.othereuk.fasta ${ouf}/euk.${hc}.seed.diamond.prokaryotes.fasta > ${ouf}/euk.${hc}.seqs.fasta

		# partition, and keep only partition that contains original OG
		echo "# mclclus"
		do_homology_clusters ${ouf}/euk.${hc}.seqs.fasta ${ouf}/euk.${hc}.seqs 4 1.2 10
		clu_list=$(awk '{ print $1 }' ${ouf}/euk.${hc}.seqs.mcl_filtered.txt | sort -u -V)

		for clu in $clu_list ; do
			awk -v clu="$clu" '$1 == clu { print $2 }' ${ouf}/euk.${hc}.seqs.mcl_filtered.txt > ${ouf}/euk.${hc}.seqs.${clu}.list
			if [ $( grep -c "^arc_\|^bac_\|^eukoth_" ${ouf}/euk.${hc}.seqs.${clu}.list) -gt 0 ] && [ $( grep -v -c "^arc_\|^bac_\|^eukoth_" ${ouf}/euk.${hc}.seqs.${clu}.list) -gt 0 ] ; then
				xargs faidx ${ouf}/euk.${hc}.seqs.fasta < ${ouf}/euk.${hc}.seqs.${clu}.list >  ${ouf}/euk.${hc}.seqs.${clu}.fasta
				echo "# launch ${ouf}/euk.${hc}.seqs.${clu}.fasta | proeuk.${hc}.${clu} | $(grep -c '>' ${ouf}/euk.${hc}.seqs.${clu}.fasta) hits ( $(grep -c '>arc_' ${ouf}/euk.${hc}.seqs.${clu}.fasta) arc, $(grep -c '>bac_' ${ouf}/euk.${hc}.seqs.${clu}.fasta) bac, $(grep -c '>eukoth_' ${ouf}/euk.${hc}.seqs.${clu}.fasta) other euk )"
				qsub -N proeuk.${hc}.${clu} -pe smp 4 $qsub_alignment ${ouf}/euk.${hc}.seqs.${clu}.fasta 4 2> /dev/null 1> /dev/null
			fi
		done

		#echo "# launch proeuk.${hc} | $(grep -c \">\" ${ouf}/euk.${hc}.seqs.fasta)"
		#qsub -N proeuk.${hc} -pe smp 4 $qsub_alignment ${ouf}/euk.${hc}.seqs.fasta 4 2> /dev/null 1> /dev/null

	done

done
rm ${ouf}/*.clstr

echo "Done"
