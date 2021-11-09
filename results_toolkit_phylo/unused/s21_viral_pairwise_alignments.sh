o="/users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/"
mkdir -p results_viruses/pairwise_alignments

for f in AF9 BIR Bromodomain Chromo CupinJmjC DOT1 GNAT_acetyltr_2 HIF-1 Hist_deacetyl Histone Kelch LinkerHistone PHD PTIP SAM SET SIR2 SNF2_N TUDOR WD40 zf-C2H2 zf-CCHH zf-CXXC Acetyltransf_1 ; do

        # internal viruses?
        # diamond blastp -q ${o}/vir.${f}.domains.fasta -d ${o}/vir.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.vir-vir.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mism$

        # all relative to euks
        echo "# align $f to euks"
        diamond blastp -q ${o}/vir.${f}.domains.fasta -d ${o}/euk.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.vir-euk.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/arc.${f}.domains.fasta -d ${o}/euk.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.arc-euk.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/bac.${f}.domains.fasta -d ${o}/euk.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.bac-euk.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos

        # all relative to bacs
        echo "# align $f to bacs"
        diamond blastp -q ${o}/vir.${f}.domains.fasta -d ${o}/bac.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.vir-bac.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/arc.${f}.domains.fasta -d ${o}/bac.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.arc-bac.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/euk.${f}.domains.fasta -d ${o}/bac.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.euk-bac.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos

        # all relative to arcs
        echo "# align $f to arcs"
        diamond blastp -q ${o}/vir.${f}.domains.fasta -d ${o}/arc.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.vir-arc.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/bac.${f}.domains.fasta -d ${o}/arc.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.bac-arc.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/euk.${f}.domains.fasta -d ${o}/arc.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.euk-arc.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos

        # all relative to virs
        echo "# align $f to virs"
        diamond blastp -q ${o}/arc.${f}.domains.fasta -d ${o}/vir.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.arc-vir.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/bac.${f}.domains.fasta -d ${o}/vir.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.bac-vir.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos
        diamond blastp -q ${o}/euk.${f}.domains.fasta -d ${o}/vir.${f}.domains.fasta -o results_viruses/pairwise_alignments/diamond.${f}.euk-vir.csv --more-sensitive --quiet -k 1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp ppos

        gzip results_viruses/pairwise_alignments/diamond.${f}*.csv

done

