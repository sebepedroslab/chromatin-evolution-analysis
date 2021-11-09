# Identification of N-terminal tails in archaeal histones

## Work

Steps:

1. Get archaeal histones from 0th connected component (HMfB-like) and create summary plots.

```bash
Rscript s01_archaeal_tails.R
```

2. Find clusters of similar tails with *MCL* (doesn't work very well for most tails, but it's a start).

```bash
makeblastdb -dbtype prot -parse_seqids -in tails_archaea.fasta
blastp -task blastp-short -query tails_archaea.fasta -db tails_archaea.fasta -out tails_to_tails.blastpshort.csv -outfmt 6 -evalue 20000
awk '$12 > 20' tails_to_tails.blastpshort.csv > tails_to_tails.blastpshort_pairs.csv
awk '{ print $1,$2,$12 }' tails_to_tails.blastpshort_pairs.csv > tails_to_tails.blastpshort_pairs.abc
mcl tails_to_tails.blastpshort_pairs.abc --abc -I 1.6 -o tails_to_tails.blastpshort_pairs.mcl.csv 2> /dev/null
```

3. Align clusters with MSA:

```bash
# Manual cluster: MXKK sequences
bioawk -c fastx '{ print $1,$4,$2 }' tails_archaea.fasta | grep "\bM.KK" | grep -v "RLI05273.1\|RLI96196.1" | awk 'BEGIN { FS="\t"  } { print ">"$1"|"$2"\n"$3 }' | tr ' ' '_' > tails_archaea-sub-MXKK.fasta
mafft --globalpair --thread 1 --reorder --maxiterate 10000 tails_archaea-sub-MXKK.fasta > tails_archaea-sub-MXKK.g.fasta

# MCL clusters: other similar pairs
bioawk -c fastx '{ print $1,$4,$2 }' tails_archaea.fasta | fgrep -f <(awk 'NR == 1' tails_to_tails.blastpshort_pairs.mcl.csv | tr '\t' '\n') | awk 'BEGIN { FS="\t"  } { print ">"$1"|"$2"\n"$3 }' | tr ' ' '_' > tails_archaea-sub-mcl1.fasta
bioawk -c fastx '{ print $1,$4,$2 }' tails_archaea.fasta | fgrep -f <(awk 'NR == 2' tails_to_tails.blastpshort_pairs.mcl.csv | tr '\t' '\n') | awk 'BEGIN { FS="\t"  } { print ">"$1"|"$2"\n"$3 }' | tr ' ' '_' > tails_archaea-sub-mcl2.fasta
bioawk -c fastx '{ print $1,$4,$2 }' tails_archaea.fasta | fgrep -f <(awk 'NR == 4' tails_to_tails.blastpshort_pairs.mcl.csv | tr '\t' '\n') | awk 'BEGIN { FS="\t"  } { print ">"$1"|"$2"\n"$3 }' | tr ' ' '_' > tails_archaea-sub-mcl4.fasta
# align
mafft --globalpair --thread 1 --reorder --maxiterate 10000 tails_archaea-sub-mcl1.fasta > tails_archaea-sub-mcl1.g.fasta
mafft --globalpair --thread 1 --reorder --maxiterate 10000 tails_archaea-sub-mcl2.fasta > tails_archaea-sub-mcl2.g.fasta
mafft --globalpair --thread 1 --reorder --maxiterate 10000 tails_archaea-sub-mcl4.fasta > tails_archaea-sub-mcl4.g.fasta

# Manual cluster: high K-frequency tails from Lokiarchaeota and Heimdallarchaeota:
head -n 26 summary_archaeal_tails.csv | awk 'NR>1' | cut -f1 > tails_archaea-topKfreq.list
xargs faidx -d ' ' tails_archaea.fasta < tails_archaea-topKfreq.list >  tails_archaea-topKfreq.fasta
mafft --globalpair --thread 1 --reorder --maxiterate 10000 tails_archaea-topKfreq.fasta > tails_archaea-topKfreq.g.fasta
```

## Data

Good hits found by [Henneman et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007582):

```bash
>arc_OLS22328.1_21-91 Archaeal histone A [Candidatus Heimdallarchaeota archaeon LC_3]
SRSKGQAFASAKVEKLIREAGAFRVSSGAIKALNDLLGERGLEVARYSVEIARNSGRRTI
KETDVALSSSK

>arc_PIN66802.1_27-91 histone [Candidatus Huberiarchaeum crystalense]
LIIPESVAIRLFKVAGAPRVSKEARDALLNLIAKYGRDVAERAVKFSKHAKRQTITSEDI
KLALE

>arc_KYH36356.1_28-105 transcription factor CBF/NF-Y/histone domain-containing protein [Candidatus Bathyarchaeota archaeon B23]
LKGSGLAEEFTLAPMRRLLKRFGELRVSVEASEELRRAVGEYGERIARAAVAHALREGRR
TVLARDVKAAREEVEGGR
```
