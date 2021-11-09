# Histone fusions

1. Find candidates

```bash
# get list of domains
cat tables/euk.Histone.domains.csv <(sed "s/^/arc_/" tables/arc.Histone.domains.csv) <(sed "s/^/bac_/" tables/bac.Histone.domains.csv) <(sed "s/^/vir_/" tables/vir.Histone.domains.csv) > histone_domain_coordinates.csv

# get domain lengths
cat <(cat tables/arc.Histone.seqs.pfamscan.csv | grep -v "^#"| grep -v "^$" | sed "s/^/arc_/") <(cat tables/euk.Histone.seqs.pfamscan.csv |  grep -v "^#"| grep -v "^$") <(cat tables/bac.Histone.seqs.pfamscan.csv | grep -v "^#"| grep -v "^$" | sed "s/^/bac_/") <(cat tables/vir.Histone.seqs.pfamscan.csv | grep -v "^#"| grep -v "^$" | sed "s/^/vir_/")| tr -s ' ' '\t' | cut -f 7,11 | sort -u > domain_lengths.csv

#cat tables/euk.Histone.domains.fasta  tables/bac.Histone.domains.fasta tables/arc.Histone.domains.fasta tables/vir.Histone.domains.fasta | bioawk -c fastx '{ print $1,$1 }' | sed -E "s/_([0-9]+)-([0-9]+)$/\t\1\t\2/" > histone_domains.csv
```

2. Annotate histone types and taxonomy for each element of the dimer:

```bash
Rscript s01_analyse_fusions.R
```

3. Establish histone fusion contiguity for species where we have genome assemblies and transcriptomes

```bash
# species list: Azfi Ddis Dpul Drer Hetalb Klenit Mgut Nvec Perkma Ppat Sphfal Symmic Ttra Apla Lgig Spur Adig Exapal Skow Ctel Spis
bash s02_txvalidation_expr-assembly_2020-12-16.sh # requires having downloaded SRA data
```

4. Annotate evidence:

```bash
cat data/fusions*.fasta > fusions_all.fasta
Rscript s03_summarise_euk_histone_fusions.R
```

5. Align consistent fusions:

```bash
Rscript s04_plot_alignments_conserved_fusions.R
```

## Old

1. Find candidates:

```bash
cut -f1,3 ../histone_classification/euk.Histone.to_histdb.Spring-Naive.csv  | sed -E "s/_[0-9]+-[0-9]+\t/\t/" | sort -k1,1| cut -f1 | uniq -c| awk '$1>1 { print $2}'| fgrep -f -  ../histone_classification/euk.Histone.to_histdb.Spring-Naive.csv | sort -k1,1 > fusion_candidates.tsv

# retrieve domains
cut -f1 fusion_candidates.tsv > fusion_candidates.domains.txt
xargs faidx -d ' ' ../histone_classification/euk.Histone.domains.fasta < fusion_candidates.domains.txt > fusion_candidates.domains.fasta

# retrieve full genes
cut -f1 fusion_candidates.tsv | sed -E "s/_[0-9]+-[0-9]+$//" | sort -u > fusion_candidates.genes.txt
xargs faidx -d ' ' ../histone_classification/euk.Histone.seqs.fasta < fusion_candidates.genes.txt > fusion_candidates.genes.fasta
```

2. How to discriminate broken from fusions?

* Easy... size!