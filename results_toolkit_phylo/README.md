# hPTM read-out toolkit

Run these commands from the present folder.

## Pfam database

Download Pfam Release 33:

```bash
mkdir -p data/databases_pfam
cd data/databases_pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.dat.gz
gunzip Pfam-A.hmm.gz
gunzip Pfam-A.hmm.dat.gz
hmmpress Pfam-A.hmm
cd ../../
```

## Taxon sampling

### Eukaryota

Ad-hoc database:

```bash
/users/asebe/xgraubove/histonome-ops/data/sequences/seq_Eukaryota.fasta
```

### Bacteria, Archaea, Viruses

1. Download protein databases from NCBI, using `esearch` and `efetch` from the NCBI Entrez Direct UNIX e-utilities.

```bash
# archaea
esearch -db protein -query "txid2157[Organism]" -sort "Taxonomy ID" | efetch -format fasta > data/databases_seqs/seq_Archaea.fasta

# virus
# add Nucleo-Cytoplasmic Large DNA Viruses from https://figshare.com/projects/NCLDV/71138
esearch -db protein -query "txid10239[Organism]" -sort "Taxonomy ID" | efetch -format fasta > data/databases_seqs/seq_Viruses.fasta

# bacteria: too many sequences in main database, so we'll have to use Refseq (release 99)
# esearch -db protein -query "txid2[Organism] AND Refseq[Filter]" -sort "Taxonomy ID" | efetch -format fasta > data/databases_seqs/seq_Bacteria.fasta
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant*.protein.faa.gz
zcat <path>/bacteria.nonredundant*.protein.faa.gz > data/databases_seqs/seq_Bacteria.fasta
```

2. Download taxonomy:

```bash
# store taxonomy (large files, ~1gb required, can't be uploaded to repo)
mkdir -p data/taxonomy/
# NCBI taxonomy from FTP server
wget -P data/taxonomy/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -xzvf data/taxonomy/new_taxdump.tar.gz -C data/taxonomy/
wget -P data/taxonomy/ https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt

# prepare taxonomy db
# field names:
# tax_id, tax_name, species, genus, family, order, class, phylum, kingdom, superkingdom
sed "s/\t|//g" data/taxonomy/rankedlineage.dmp > data/taxonomy_ranked.tsv && gzip data/taxonomy_ranked.tsv


# species codes to txid
# wget -P tmp https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# tar -xzvf tmp/taxdump.tar.gz -C tmp/
# txid to TOL domain (archaea, eukaryota, virus...)
# wget -P tmp https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz
# tar -xzvf tmp/taxcat.tar.gz -C tmp/
```

3. Format database and taxonomy:

```bash
# extract taxonomic info from FASTA, using exact species names
tax="Bacteria" # same with Archaea, Viruses and Bacteria
db=data/databases_seqs/seq_${tax}.fasta
paste <(bioawk -c fastx '{ print $1}' ${db}) <(bioawk -c fastx '{ print $4}' ${db} | awk -F  "[" '{ if (NF>1) { print $NF } else { print "" }}' | sed "s/]//" ) > ${db}.taxa.csv
gzip ${db}.taxa.csv

# also, prepare blast db (used for sequence extraction)
makeblastdb -dbtype prot -parse_seqids -in ${db}
```

## Gene family analyses

Format input databases:

```bash
esl-sfetch --index seq_Eukaryota.fasta
```

HMM searches of each profile in each database (eukaryotes, bacteria, archaea, viruses...).

```bash
# list of families, HMMs and inflation values used in each search are available here:
# data/gene_families_hmm.csv
bash s01_hmmsearches_v07_2020-06-23.sh <gene family> <comma-separated list of hmms> <inflation> <genomes fasta> <do phylogenies? Y/n> <prefix> <output folder alignments> <output folder hmmsearches> GA
# calls `qsub_alignment.sh` and `qsub_pfamscan.sh` scripts, that contain the specific 
# commands used to compute alignments, trees etc.


# specifically, for each taxon set:
# euk
while read -a i ; do bash s01_hmmsearches_v08_2020-11-12.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ~/histonome-ops/data/sequences/seq_Eukaryota.fasta Y euk /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/ ${i[4]} ; done < ../data/gene_families_hmm.csv
# vir
while read -a i ; do bash s01_hmmsearches_v08_2020-11-12.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ~/histonome-ops/data/sequences/seq_Virus.fasta N vir /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/ ${i[4]} ; done < ../data/gene_families_hmm.csv
# arc
while read -a i ; do bash s01_hmmsearches_v08_2020-11-12.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ~/histonome-ops/data/sequences/seq_Archaea.fasta N arc /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/ ${i[4]} ; done < ../data/gene_families_hmm.csv
# bac
while read -a i ; do bash s01_hmmsearches_v08_2020-11-12.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ~/histonome-ops/data/sequences/seq_Bacteria.fasta N bac /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/ ${i[4]} ; done < ../data/gene_families_hmm.csv

# for bacterial, archaeal and viral histones, launch pfamscan searches:
qsub -N pf.virhis qsub_pfamscan.sh ../../results-phylogenies/searches_Nov20/vir.Histone.seqs.fasta  /users/asebe/xgraubove/data/pfam/ 4
qsub -N pf.bachis qsub_pfamscan.sh ../../results-phylogenies/searches_Nov20/bac.Histone.seqs.fasta  /users/asebe/xgraubove/data/pfam/ 4
qsub -N pf.archis qsub_pfamscan.sh ../../results-phylogenies/searches_Nov20/arc.Histone.seqs.fasta  /users/asebe/xgraubove/data/pfam/ 4
```

**WARNING**:

* Kill jobs for WD40, zf-C2H2 and Kelch, and launch dedicated diamond searches using human seeds for selected proteins instead:

```bash
diamond blastp -d ~/histonome-ops/data/sequences/seq_Eukaryota.fasta  -q seed_fasta/list_Kelch_human.fasta -o list_Kelch_human.diamond.csv --more-sensitive  -k 1000 --quiet
diamond blastp -d ~/histonome-ops/data/sequences/seq_Eukaryota.fasta  -q seed_fasta/list_WD40_human.fasta -o list_WD40_human.diamond.csv --more-sensitive  -k 1000 --quiet
diamond blastp -d ~/histonome-ops/data/sequences/seq_Eukaryota.fasta  -q seed_fasta/list_C2H2_human.fasta -o list_C2H2_human.diamond.csv --more-sensitive  -k 1000 --quiet
diamond blastp -d ~/histonome-ops/data/sequences/seq_Eukaryota.fasta  -q seed_fasta/list_PCRing.fasta -o list_PCRing.diamond.csv --more-sensitive  -k 1000 --quiet
# separate by seed protein and launch independent trees with qsub
# while read i ; do awk '$11 < 1e-10' list_Kelch_human.diamond.csv | grep -w $i  | cut -f2 | sort -u > list_Kelch_seqs.$i.txt ; done < list_Kelch_human.txt
# for RINGs:
awk '$11 < 1e-10' seed_fasta/list_PCRing.diamond.csv | cut -f2 | sort | uniq -c | sort -n| awk '$1 > 1 { print $2 }'| sort > seed_fasta/list_PCRing.seqs.txt
xargs faidx -d ' ' ../../data/sequences/seq_Eukaryota.fasta < seed_fasta/list_PCRing.seqs.txt > seed_fasta/list_PCRing.seqs.fasta
# now copy these output files
```

Once all this work is done, pull the necessary data into the repository:

```bash
# add data to gene_counts/ folder:
for s in euk arc bac vir ; do while read -a i ; do cat ${s}.${i[2]}.genes.list | sed "s/$/\t${i[2]}/" ; done < ../../histonome-analysis/data/gene_families_hmm.csv  > /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_toolkit_phylo/gene_counts/${s}_genecounts.csv ; done

# add data to gene_trees/
cp /nfs/users2/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/*treefile /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_toolkit_phylo/gene_trees

# add gene sequences to gene_sequences/
cp /nfs/users2/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_Nov20/*.seqs.fasta /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_toolkit_phylo/gene_sequences

# add architectures to architectures.csv
cat /nfs/users2/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/euk.*.seqs.pfamscan_archs.csv | sort -u > /users/asebe/xgraubove/histonome-ops/histonome-analysis/results_toolkit_phylo/architectures.csv
```

For the eukaryotic dataset, call orthologous groups from phylogenies using species overlap (results in `results_evol`):

```bash
# call orthologs from each tree:
for i in gene_trees/*treefile ; do python s02_parse_phylogeny_2020-11-30.py -i $i -o gene_trees/ -r human_gene_names.csv  -refsps Hsap -itermidroot 10 -min_transfer_support 50 -cut_gene_names 100 -ogprefix $(basename $i | cut -f2,3 -d '.'). -p $(basename $i | sed "s/.seqs.iqtree.treefile/.possom/") ; done

# concatenate all OGs in a single table:
# beware: step to remove Plebac, which is completely synonymous with Pbla in trees
echo -e "gene\torthogroup\torthologous_to" > orthogroups_euk.csv
for i in gene_trees/*ortholog_groups.csv ; do awk 'NR > 1' $i | grep -v "^Plebac_" ; done >> orthogroups_euk.csv
```

## Ancestral reconstruction

Ancestral reconstruction with **Count**:

1. Sequential training:

```bash
# get training dataset (from variable PFAM families)
Rscript s03_ancestral_reconstruction-prepare.R

# # one gamma cat
# java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -gain_k 1 -loss_k 1 -uniform_duplication true  data/species_tree.newick data/pfam_domain_counts.train.csv > data/pfam_domain_counts.rates.r100.g1.txt
# # two gamma cats
# java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -gain_k 2 -loss_k 2 -uniform_duplication true  data/species_tree.newick data/pfam_domain_counts.train.csv data/pfam_domain_counts.rates.r100.txt > data/pfam_domain_counts.rates.r100.g2.txt

### gamma categories for gain, loss, length and transfer (transfer seems important because it emulates "errors")
### THIS IS SLOW AS HELL, RERUN WITH CARE
# start without gamma cats
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -uniform_duplication true  data/species_tree.newick data/pfam_domain_counts.train.csv > data/pfam_domain_counts.rates.G0.txt
# add one gamma cat everywhere
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -gain_k 1 -loss_k 1 -transfer_k 1 -length_k 1 -uniform_duplication true  data/species_tree.newick data/pfam_domain_counts.train.csv data/pfam_domain_counts.rates.G0.txt > data/pfam_domain_counts.rates.G1.txt
# two gamma cats
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.ML -opt_rounds 100 -gain_k 2 -loss_k 2 -transfer_k 2 -length_k 2 -uniform_duplication true  data/species_tree.newick data/pfam_domain_counts.train.csv data/pfam_domain_counts.rates.G1.txt > data/pfam_domain_counts.rates.G2.txt
```

2. Run posterior probability analysis:

```bash
java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 data/species_tree.newick orthogroups_euk.ancestral.counts.csv data/pfam_domain_counts.rates.G2.txt  > orthogroups_euk.ancestral.posteriors.csv

# based on presence only? This would require training based on pres/abs, not counts, so NO
# java -cp ~/Programes/Count/Count.jar  ca.umontreal.iro.evolution.genecontent.Posteriors -max_paralogs 10000 data/species_tree.newick orthogroups_euk.ancestral.presence.csv data/pfam_domain_counts.rates.G2.txt  > orthogroups_euk.ancestral.presence_posteriors.csv

# in addition, run Wagner parsimony with various gain/loss ratios (classical Wagner)
# If g>1, scattered distributions are explained by multiple losses (more presence at LECA)
# If g<1, there'll be multiple gains (i.e., lateral transfers, less presence at LECA)
java -cp ~/Programes/Count/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner -gain 5 -max_paralogs 10000 data/species_tree.newick orthogroups_euk.ancestral.counts.csv > orthogroups_euk.ancestral.wagner_g5.csv
# grep FAMILY  presence lines
grep "# FAMILY" orthogroups_euk.ancestral.wagner_g5.csv > orthogroups_euk.ancestral.wagner_g5_pres.csv

```

3. Ancestral reconstructions:

```bash
# plots of gene family evolutionary history, and heatmaps:
Rscript s04a_ancestral_reconstruction-posterior.R

# pie plots with evidence of presence of OGs at LECA:
Rscript s04b_LECA_evidence_2021-01-28.R

# pie plots with evidence of presence of OGs gained at the last Metazoan common ancestor:
Rscript s04c_LMetCA_evidence_20201-09-16.R
```

4. Create annotated tables of orthologs, including data for reference species:

```bash
# get reference species
cd data/
i="Scer"
f=ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz

i="Atha"
f=ftp://ftp.ensemblgenomes.org/pub/release-47/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz

i="Dmel"
f=ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz

# do
fa=$(basename ${f}) ;  wget ${f} ; gunzip ${fa} ; paste <(grep "gene_symbol" ${fa%%.gz} | awk '{ print $1 }'| sed "s/>/${i}_/") <(grep -o "gene_symbol:[^ ]*" ${fa%%.gz} | sed "s/gene_symbol://") > gene_names_${i}.csv ; rm ${fa%%.gz}


```

### Old commands

Previous pipeline, based on Dollo parsimony. Scripts in `unused`:

```bash
# # launch species tree-based evolutionary reconstruction
# # Rscript s03a_orthogroup_consolidation.R # this fuses orthogroups that are "islike" a single master OG
# python s03b_possom_reconstruction.py -tree data/species_tree.newick -ort orthogroups_euk.csv -out orthogroups_euk.possom

# # plots of gene family evolutionary history, and heatmaps:
# Rscript s04a_orthogroup_along_phylogeny_2020-10-29.R

# # pie plots with evidence of presence of OGs at LECA:
# s04b_LECA_evidence_2021-01-28.R
```

## Gene family counts

For all datasets, count number of hits for each gene family (results in `results_counts`):

```bash
# for eukaryotes
Rscript s05a_gene_counts_hmm-euk_2020-10-29.R

# for non-eukaryotes:
Rscript s05b_gene_counts_hmm-noneuk_2020-10-29.R
# WARNING:
# bacteria NCBI taxonomy is huge: subset it to the genes present in our searches (it doesn't matter if other species go missing, because for counting purposes these are pulled from the general taoxnomy file)
# fgrep -w -f <(cut -f1 gene_counts/bac_genecounts.csv ) <(zcat /users/asebe/xgraubove/histonome-ops/data/sequences/seq_Bacteria.taxa.csv.gz) > ../data/seq_Bacteria.taxa.csv && gzip ../data/seq_Bacteria.taxa.csv
# alternatively: copy it manually from cluster! can't be included in the repository
# HERE: /users/asebe/xgraubove/histonome-ops/data/sequences/seq_Bacteria.taxa.csv.gz

# correlation of histone presence with other enzymes, in non-euks
s05c_histone_correlation_noneuk_2021-03-25.R
```

For the eukaryotic datasets, analyse protein domain architectures in each orthogroup (results in `results_domains`):

```bash
# plot gene architecture networks
python s06_architecture_networks_2020-11-30.py
# plot architecture network for ALL READERS at the same time
python s06_architecture_networks-READERS_2021-03-15.py

# plot heatmaps of common architectures per family, OG, etc
Rscript s07_architecture_heatmaps-cooc.R
# identify matches between MOST COMMON architectures and OGs
Rscript s08_architecture_consensuses.R
```

## TE fusion analyses

First, annotate TEs in the protein sequences (and map protein sequences to assemblies if possible?).

```bash
# change dir to somewhere else
# Download Dfam 3.3 (outside of repository)
wget https://www.dfam.org/releases/Dfam_3.3/families/Dfam.hmm.gz  # for headers
wget https://www.dfam.org/releases/Dfam_3.3/families/Dfam.embl.gz # for consensus seqs

# clean EMBL file and convert to fasta:
gunzip Dfam.embl.gz
grep -v "^CC" Dfam.embl > Dfam_clean.embl
seqret -sequence Dfam_clean.embl -sformat1 embl -outseq Dfam_clean.fasta
bioawk -c fastx '{ print $1,$2,$4 }' Dfam_clean.fasta | awk '{ print ">"$3"\n"$2 }' > Dfam_clean.fasta.2 && mv Dfam_clean.fasta.2 Dfam_clean.fasta
makeblastdb -dbtype nucl -parse_seqids -in Dfam_clean.fasta

# get TE annotation taxonomy
grep "^NAME\|^ACC\|^CT" <(zcat Dfam.hmm.gz) | awk '{ print $2 }'| awk 'NR%3{printf "%s ",$0;next;}1' | awk '{ print $2,$1,$3 }'  | tr ' ' '\t' > Dfam.TEclasses.txt
cp Dfam.TEclasses.txt ~/Documents/histonome-analysis/data # path to data folder in repository

# first, prepare list of genes of interest, and dedicated fasta file:
cat results_domains/*network_genes.csv | fgrep -w -f <(cut -f1 ../data/transposable_element_domains.csv) | sort -u > results_TEfusions/architectures_with_TEs.csv
cut -f1 results_TEfusions/architectures_with_TEs.csv | sort -u  > results_TEfusions/architectures_with_TEs.genes.txt

# align TEs with candidate sequences (go back to main directory)
cat gene_sequences/*.fasta | bioawk -c fastx '{ print $1,$2 }' | fgrep -w -f results_TEfusions/architectures_with_TEs.genes.txt | sort -u | awk '{ print ">"$1"\n"$2 }' > results_TEfusions//all_hits.fasta
tblastn -db /home/xavi/dades/Dfam_clean.fasta -query results_TEfusions//all_hits.fasta -out  results_TEfusions//all_hits.Dfam.tblastn.csv -outfmt 6 -num_threads 8 -max_target_seqs 10
```

Second, validate structure of fusion genes. This consists of various steps:

1. Find genes models built across stretches of `NNNN` nucleotides (in the genome assemblies)
2. Identify genes whose predicted CDS has continuous coverage between the focus domain and the TE domain
3. Identify false positive fusions and species prone to these sort of problems.

```bash
# for each sequence with TE domains, find equivalent gene ID in GTF-able db (BLAST)
# next, RNA-seq data from SRA & map it (salmon) for all species specified in the ../data/validation_species.csv file
bash s10_txvalidation_sra-map_2020-12-16.sh # this will take a while...
# next, check expression status, gene structure and assembly contiguity for all candidate genes
bash s11_txvalidation_expr-assembly_2020-12-16.sh
```

In parallel, launch TE phylogenies:

```bash
# launch phylogenies (takes a while...)
while read -a i ; do bash s01_hmmsearches_v08_2020-12-21-DOMAINS.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ~/histonome-ops/data/sequences/seq_Eukaryota.fasta Y euk /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_TEs_Dec20/ ${i[4]} ; done < ../data/transposon_search_hmm.csv

# find TE clusters
mkdir -p gene_trees_TE/
cp  /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/euk*treefile gene_trees_TE/
for i in gene_trees_TE/*treefile ; do python s02_parse_phylogeny_2020-11-30.py -i $i -o gene_trees_TE/ -r human_gene_names.csv  -refsps Hsap -itermidroot 10 -min_transfer_support 50 -cut_gene_names 100 -ogprefix $(basename $i | cut -f2,3 -d '.'). -p $(basename $i | sed "s/.seqs.iqtree.treefile/.possom/") ; done

# # dedicated domain-specific phylogenies for Chromo, SNF and PHD domains
# # THIS IS DEPRECATED
# bash s01_hmmsearches_v08_2020-12-21-DOMAINS.sh Chromo Chromo,Chromo_shadow 1.1 2 ~/histonome-ops/data/sequences/seq_Eukaryota.fasta Y cor /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_TEs_Dec20/ GA
# bash s01_hmmsearches_v08_2020-12-21-DOMAINS.sh PHD PHD,PHD_2,PHD_3,PHD_4 1.1 2 ~/histonome-ops/data/sequences/seq_Eukaryota.fasta Y cor /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_TEs_Dec20/ GA
# bash s01_hmmsearches_v08_2020-12-21-DOMAINS.sh SNF2_N SNF2_N 1.1 2 ~/histonome-ops/data/sequences/seq_Eukaryota.fasta Y cor /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_TEs_Dec20/ GA
# find OGs
# mkdir -p gene_trees_TEcore/
# cp  /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_TEs_Dec20/cor*treefile gene_trees_TEcore/
# for i in gene_trees_TEcore/cor*treefile ; do python s02_parse_phylogeny_2020-11-30.py -i $i -o gene_trees_TEcore/ -r human_gene_names.csv  -refsps Hsap -itermidroot 10 -min_transfer_support 50 -cut_gene_names 100 -ogprefix $(basename $i | cut -f2,3 -d '.'). -p $(basename $i | sed "s/.seqs.iqtree.treefile/.possom/") ; done
```

Finally, analyse candidate TE fusions:

1. Record expression status and (if available) stage-specificity.
2. Record assembly contiguity
3. Record gene structure: monoexonic points to retrocopies and possibly machinery hijacking, multiexonic might point to TE domestication
4. Record TE class (Dfam analysis above).
5. Record TE position in the genome to evaluate if expansions are due to retrocopies (scattered) or tandem duplications (clustered)
6. Add info from TE phylogenies

```bash
Rscript s12_TEannot_2021-01-19.R
```

Align `DUF1087` with Dfam hits:

```bash
# grep -w "DUF1087" results_TEfusions/gene_fusion_evidence.csv | cut -f1 > results_TEfusions/fusions_DUF1087.seqs_euks.list
grep -w "DUF1087" results_TEfusions/gene_fusion_evidence.csv | cut -f2 > results_TEfusions/fusions_DUF1087.seqs_dfam.list
cat results_TEfusions/data/fusions.*.seqs.fasta | bioawk -c fastx '$1,$2'| grep "DUF1087" | awk '{print ">"$1"\n"$2 } ' > results_TEfusions/fusions_DUF1087.seqs_euks.fasta
bioawk -c fastx '{ print $1 }' results_TEfusions/fusions_DUF1087.seqs_euks.fasta| cut -f1 -d '|' | sed "s/\d+-\d+//" | sort -u > results_TEfusions/fusions_DUF1087.seqs_euks.list
# get dfam blast hits
tblastn -db /home/xavi/dades/Dfam_clean.fasta -query results_TEfusions/fusions_DUF1087.seqs_euks.fasta -out  results_TEfusions/fusions_DUF1087.tblastn_dfam.csv -outfmt 6 -num_threads 8 -max_target_seqs 10
cut -f2 results_TEfusions/fusions_DUF1087.tblastn_dfam.csv | sort -u > results_TEfusions/fusions_DUF1087.tblastn_dfam.list
xargs faidx -d ' ' /home/xavi/dades/Dfam_clean.fasta < results_TEfusions/fusions_DUF1087.tblastn_dfam.list > results_TEfusions/fusions_DUF1087.tblastn_dfam.fasta
# get ORFs from dfam
getorf results_TEfusions/fusions_DUF1087.tblastn_dfam.fasta -outseq results_TEfusions/fusions_DUF1087.orf_dfam.fasta
sed "s/ \[/_/" results_TEfusions/fusions_DUF1087.orf_dfam.fasta | tr -d ' ]' | sed "s/(REVERSESENSE)//" > results_TEfusions/fusions_DUF1087.orf_dfam.fasta.2 && mv results_TEfusions/fusions_DUF1087.orf_dfam.fasta.2 results_TEfusions/fusions_DUF1087.orf_dfam.fasta
# reblst on ORGs
makeblastdb -dbtype prot -parse_seqids -in results_TEfusions/fusions_DUF1087.orf_dfam.fasta
blastp -db results_TEfusions/fusions_DUF1087.orf_dfam.fasta -query results_TEfusions/fusions_DUF1087.seqs_euks.fasta -out  results_TEfusions/fusions_DUF1087.orf_dfam.csv -outfmt 6 -num_threads 8 -max_target_seqs 10
cut -f2 results_TEfusions/fusions_DUF1087.orf_dfam.csv > results_TEfusions/fusions_DUF1087.orf_dfam.list
xargs faidx -d ' ' results_TEfusions/fusions_DUF1087.orf_dfam.fasta < results_TEfusions/fusions_DUF1087.orf_dfam.list > results_TEfusions/fusions_DUF1087.orf_dfam_filt.fasta

# concatenate
cat results_TEfusions/fusions_DUF1087.seqs_euks.fasta results_TEfusions/fusions_DUF1087.orf_dfam_filt.fasta > results_TEfusions/fusions_DUF1087.seqs_eukdfam.fasta

# launch tree
bash qsub_alignment-single.sh results_TEfusions/fusions_DUF1087.seqs_eukdfam.fasta 10
```

Currently missing information:

* Identify signals of purifying selection in the TE domain that might indicate domestication (`dN/dS`), as in SETMAR.

## Viral toolkit analysis

1. Concatenate all hits and launch phylogenies of homology groups with eukaryotic and viral sequences:

```bash
# launch phylogenies with viral sequences (single-run tree)
bash s21_viralsearches_2021-01-14.sh
Rscript s22_viral_annotation-phylogeny_2021-01-15.R # this is a bit slow because the Acetyltransferase phylogeny is big but hey
```

2. Virus Pfam domain architectures:

```bash
# get non-eukaryotic pfam architectures
i=vir
cat /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/vir.*.seqs.fasta | bioawk -c fastx '{ print $1,$2 }' | sort -u | awk '{ print ">'$i'_"$1"\n"$2 }' > gene_sequences/noneuk_vir.fasta
qsub -pe smp 12 -N noneuk_vir.fasta qsub_pfamscan.sh gene_sequences/noneuk_vir.fasta /users/asebe/xgraubove/data/pfam
```

3. Parse annotations of viral homologs:

```bash
Rscript s23_viral_annotation-post_2021-02-26.R
```

## Prokaryotic roots of the eukaryotic toolkit

For each eukaryotic OG, we'll retrieve 1) the OG, 2) its closest prokaryotic homologs, and 3) other eukaryotic homologs, and lunch as many trees as warranted by MCL partitioning. Then, we'll load the tree collection and ascertain whether the OG of interest is most similar to its eukaryotic counterparts (indicating intra-eukaryotic emergence of that family) or to either archaea or bacteria. If warranted, we'll record the taxonomic profile of the prokaryotic homologs. That way, we'll be able to identify OGs that have clear prokaryotic roots.
This information can be aggregated at the gene family level.

1. Obtain trees with the prokaryotic homologs of each eukaryotic homology group (only selected gene families present in archaea?).

```bash
bash s41_preeuk_search_famcentric_2021-03-21.sh
```

2. Parse OG-centric phylogenies:

```bash
Rscript s42_preeuk_annotation_famcentric-phylogeny_2021-03-22.R
```

## Archaeal toolkit analysis

WARNING: This may end up being useless.

1. Concatenate all hits and launch phylogenies of homology groups with eukaryotic and archaeal sequences:

```bash
# launch phylogenies with archaeal sequences (single-run tree)
bash s31_archaeal_search_2021-03-18.sh
Rscript s32_archaeal_annotation-phylogeny_2021-03-18.R
```

2. Archaea Pfam domain architectures:

```bash
# get non-eukaryotic pfam architectures
cat /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_Nov20/arc.*.seqs.fasta | bioawk -c fastx '{ print $1,$2 }' | sort -u | awk '{ print ">arc_"$1"\n"$2 }' > gene_sequences/noneuk_arc.fasta
qsub -pe smp 12 -N noneuk_arc.fasta qsub_pfamscan.sh gene_sequences/noneuk_arc.fasta /users/asebe/xgraubove/data/pfam
```

3. Parse annotations of archaeal homologs:

```bash

### ADAPT:
Rscript s33_archaeal_annotation-post_2021-03-18.R
```

The end.

```python
        ~+ the end +~

                 *       +
           '                  |
       ()    .-.,="``"=.    - o -
             '=/_       \     |
          *   |  '=._    |
               \     `=./`,        '
            .   '=.__.=' `='      *
   +                         +
        O      *        '       .
```

## New?

Joint analysis of histones and ribosomal proteins?

```bash
# all together
bash s01_hmmsearches_v08_2020-11-12.sh histplus Ribosomal_S6e,Histone,CBFD_NFYB_HMF,CENP-T_C 1.1 10  ~/histonome-ops/data/sequences/seq_Eukaryota.fasta N euk /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_HistRib_Abr21/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_HistRib_Abr21/ 0.001
bash s01_hmmsearches_v08_2020-11-12.sh histplus Ribosomal_S6e,Histone,CBFD_NFYB_HMF,CENP-T_C 1.1 10  ~/histonome-ops/data/sequences/seq_Archaea.fasta N arc /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_HistRib_Abr21/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_HistRib_Abr21/ 0.001

# ribosomal s6e
bash s01_hmmsearches_v08_2020-11-12.sh Ribosomal_S6e  Ribosomal_S6e 1.1 10  ~/histonome-ops/data/sequences/seq_Archaea.fasta N arc /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_HistRib_Abr21/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_HistRib_Abr21/ 0.001
bash s01_hmmsearches_v08_2020-11-12.sh Ribosomal_S6e Ribosomal_S6e 1.1 10  ~/histonome-ops/data/sequences/seq_Eukaryota.fasta N euk /users/asebe/xgraubove/histonome-ops/results-phylogenies/alignments_HistRib_Abr21/ /users/asebe/xgraubove/histonome-ops/results-phylogenies/searches_HistRib_Abr21/ 0.001

```
