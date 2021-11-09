#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m abe
#$ -q long-sl7
#$ -l virtual_free=50G,h_rt=2592000
#$ -o tmp/
#$ -e tmp/

esl-sfetch --index /users/asebe/xgraubove/histonome-ops/data/databases_seqs/seq_Bacteria.fasta
