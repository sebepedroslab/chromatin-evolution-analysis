#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m abe
#$ -q long-sl7
#$ -l virtual_free=50G,h_rt=2592000
#$ -o tmp/
#$ -e tmp/

hmmsearch --domtblout /users/asebe/xgraubove/histonome-ops/results_toolkit_hmm/searches/hmmsearch.Bacteria.domain_CENP-T_C.domtable \
--domE 0.001 --cpu 8 \
/users/asebe/xgraubove/histonome-ops/results_toolkit_hmm/hmms/CENP-T_C.hmm \
../data/databases_seqs//seq_Bacteria.fasta

mcl /users/asebe/xgraubove/histonome-ops/results_toolkit_hmm/alignments/fam.Bacteria.Acetyltransf_1.diamond.abc \
--abc -I 1.1 \
-o /users/asebe/xgraubove/histonome-ops/results_toolkit_hmm/alignments/fam.Bacteria.Acetyltransf_1.mcl.csv
