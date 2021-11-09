# input variables
pfam_db="/users/asebe/xgraubove/data/pfam" # where to find PFAM-A.hmm database (preformatted)
qsub_alignment="qsub_alignment.sh" # qsub script to submit alignment & tree jobs
qsub_pfamscan="qsub_pfamscan.sh" # qsub script to submit pfamsan jobs
n_cpu="8"

if [ -z ${10} ] ; then
	echo "missing arguments:"
	echo "1.  family id, e.g. homeobox"
	echo "2.  comma-sparated list of HMM files, e.g. Homeodomain,Homeobox_KN"
	echo "3.  inflation for MCL, e.g. 1.1"
	echo "4.  min alignment size for phylogenetic analysis, e.g. 10"
	echo "5.  database to search, e.g. data/genomes.fasta"
	echo "6.  Do MCLclust+MSA+phylogenies? Y/N"
	echo "7.  prefix"
	echo "8.  alignments folder"
	echo "9.  searches folder"
	echo "10. evalue threshold or GA for hmmsearch"
	exit
fi

# family ids
fam_id=$1

# list of hmms belonging to this gene family
hmm_syn=$(echo $2 | tr ',' ' ')

# search parameters
inflation=$3
min_alignment_size=$4
dbfile=$5

# do phylogenies+MSA?
dophy=$6

# prefix for search files
prefix=$7

# mkdirs
alignments=$8
searches=$9
mkdir -p ${alignments}
mkdir -p ${searches}

# functions
function do_hmmsearch {

    # input
    local hmm=$1
    local out=$2
    local fas=$3
    local cpu=$4
    local thr=$5

    if [ $thr == "GA" ] ; then
    # search with gathering threshold
    hmmsearch \
        --domtblout ${out}.domtable \
        --cut_ga \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    else
    # search with domain evalue
    hmmsearch \
        --domtblout ${out}.domtable \
        --domE ${thr} \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    fi

    # clean domtable
    grep -v "^#" ${out}.domtable \
    | awk 'BEGIN { OFS="\t" } { print $1,$18,$19,$4,$5,$12 }' \
    > ${out}.domtable.csv

}

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

### RUN ANALYSES

    echo "# ${dbfile}: $fam_id | HMM search" | tee -a ${0}.${prefix}.log
    for hmm in ${hmm_syn} ; do
        do_hmmsearch hmms/${hmm}.hmm ${searches}/${prefix}.hmmsearch.domain_${hmm} ${dbfile} ${n_cpu} GA
        cat ${searches}/${prefix}.hmmsearch.domain_${hmm}.domtable.csv \
        >> ${searches}/${prefix}.${fam_id}.domains.csv.tmp
    done

    # any hits?
    cut -f1 ${searches}/${prefix}.${fam_id}.domains.csv.tmp | sort -u > ${searches}/${prefix}.${fam_id}.genes.list
    num_hit_genes=$(cat ${searches}/${prefix}.${fam_id}.genes.list | wc -l)
    echo "# ${dbfile}: ${fam_id} | # genes found = $num_hit_genes" | tee -a ${0}.${prefix}.log

    if [ ${num_hit_genes} -eq 0 ] ; then

        echo "# ${dbfile}: ${fam_id} | Omit downstream analyses" | tee -a ${0}.${prefix}.log

    else

        # merge overlapping domains
        bedtools merge \
        -i <(sort -k1,1 -k 2,2n ${searches}/${prefix}.${fam_id}.domains.csv.tmp) \
        -c 4 -o collapse \
        -d 0 \
        > ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 \
        && rm -f ${searches}/${prefix}.${fam_id}.domains.csv.tmp

        # report
        num_unique_domns=$(cut -f1 ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 | wc -l)
        echo "# ${dbfile}: ${fam_id} | # unique domains = $num_unique_domns" | tee -a ${0}.${prefix}.log

        # extract complete sequences
        esl-sfetch -f ${dbfile} ${searches}/${prefix}.${fam_id}.genes.list \
        > ${searches}/${prefix}.${fam_id}.seqs.fasta 2>> ${0}.${prefix}.log

        # dict sequence lengths
        samtools faidx ${searches}/${prefix}.${fam_id}.seqs.fasta

        # expand domain region by a fixed amount of aa
        bedtools slop \
            -i ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 \
            -g ${searches}/${prefix}.${fam_id}.seqs.fasta.fai \
            -b 10 \
            | awk '{ print $1, $2+1, $3, $4 }' \
            > ${searches}/${prefix}.${fam_id}.domains.csv \
            && rm ${searches}/${prefix}.${fam_id}.domains.csv.tmp2

        # extract domain region
        esl-sfetch --index ${searches}/${prefix}.${fam_id}.seqs.fasta 1>> ${0}.${prefix}.log 2>> ${0}.${prefix}.log
        esl-sfetch -C -f   ${searches}/${prefix}.${fam_id}.seqs.fasta \
        <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' ${searches}/${prefix}.${fam_id}.domains.csv) \
        > ${searches}/${prefix}.${fam_id}.domains.fasta 2>> ${0}.${prefix}.log

        if [ $dophy == "Y" ] ; then

        # find clusters of homology (using whole genes or domains?)
        echo "# ${dbfile}: ${fam_id} | pairwise alignments + MCL..." | tee -a ${0}.${prefix}.log
        do_homology_clusters \
            ${searches}/${prefix}.${fam_id}.domains.fasta \
            ${alignments}/${prefix}.${fam_id} \
            ${n_cpu} \
            ${inflation} \
            ${min_alignment_size} 2>> ${0}.${prefix}.log

        # report
        num_homgroups=$( awk 'END { print NR }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.csv)
        num_seqs_in_hg=$(awk 'END { print NR }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt)
        echo "# ${dbfile}: ${fam_id} | homology groups = $num_homgroups ($num_seqs_in_hg / $num_hit_genes)" | tee -a ${0}.${prefix}.log

        # pfamscan
        jobname=pf.${fam_id}
        echo "# ${dbfile}: ${fam_id} | Submit pfamscan job | $jobname" | tee -a ${0}.${prefix}.log
        qsub -N ${jobname} $qsub_pfamscan ${searches}/${prefix}.${fam_id}.seqs.fasta ${pfam_db} ${n_cpu} 1>> ${0}.${prefix}.log 2>> ${0}.${prefix}.log

        # obtain lists of clusters
        clu_list=$(awk '{ print $1 }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt | sort -u -V)

        for clu in $clu_list; do

            # retrieve sequence list
            awk -v clu="$clu" '$1 == clu { print $2 }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt \
            > ${alignments}/${prefix}.${fam_id}.${clu}.seqs.list

            # retrieve sequences from original fasta
            xargs faidx ${searches}/${prefix}.${fam_id}.domains.fasta \
            < ${alignments}/${prefix}.${fam_id}.${clu}.seqs.list \
            > ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta \
            2>> ${0}.${prefix}.log

            # report
            num_seqs_in_ali=$(grep -c ">" ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta)
            echo "# ${dbfile}: ${fam_id} | $clu | create fasta n = $num_seqs_in_ali" | tee -a ${0}.${prefix}.log

            # submit alignment job
            jobname=${fam_id}.${clu}
            echo "# ${dbfile}: ${fam_id} | $clu | submit alignment job $jobname" | tee -a ${0}.${prefix}.log
            qsub -N ${jobname} -pe smp ${n_cpu} $qsub_alignment ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta ${n_cpu} 1>> ${0}.${prefix}.log 2>> ${0}.${prefix}.log
        done # end FOR loop that iterates over homology groups


    fi # end IF statement that allows you to skip phylogenies
    fi # end IF statement that excludes empty searches


echo "# END" | tee -a ${0}.${prefix}.log

exit 0
