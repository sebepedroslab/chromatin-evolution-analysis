# Ancestral reconstruction

Reconstruct LECA using a [Scrollsaw-like](https://www.biorxiv.org/content/10.1101/2020.01.27.920793v2) approach.

Steps:

1. Reduce datasets (takes classified sequences from `../results_toolkit_phylo/gene_sequences/`):

```bash
mkdir -p reduced_datasets
Rscript s01_reduce_scrollsaw_2021-02-01.R
```

2. Launch phylogenies:

```bash
mkdir -p reduced_alignments
cp reduced_datasets/*.fasta reduced_alignments/
for i in alignments_reduced/reduced.*genes.fasta ; do qsub -N $(basename $i) -l h_rt=432000,h_vmem=60G,virtual_free=60G -pe smp 6 ../results_toolkit_phylo/qsub_alignment.sh $i 6 ; done
# wait...
cp reduced_alignments/*.treefile reduced_trees/
```

3. Get orthogroups (using taxonomic groupings instead of species in `possom`) and annotate them using pre-computed OGs (from `results_toolkit_phylo`:)

```bash
# create aOG tables (ancestral orthogroups):
for i in reduced_trees/reduced.*.treefile ; do python ../results_toolkit_phylo/s02_parse_phylogeny_2020-11-30.py -i $i -o reduced_trees/ -itermidroot 10 -split "|" -ogprefix aOG.$(basename $i | cut -f2 -d '.') -p $(basename $i | sed "s/.genes.iqtree.treefile/.possom/") ; done
# without iterative midrooting?
for i in reduced_trees/reduced.*.treefile ; do python ../results_toolkit_phylo/s02_parse_phylogeny_2020-11-30.py -i $i -o reduced_trees/ -split "|" -ogprefix aOG.$(basename $i | cut -f2 -d '.') -p $(basename $i | sed "s/.genes.iqtree.treefile/.possom/") ; done

# annotate aOGs to OGs from the full-sampling phylogenetic analysis:
Rscript s03_annotate_aOGs_2021-02-01.R
```

The end:

```perl
                                  \
                                  `\,/
                                  .-'-.
                                 '     `
                                 `.   .'
                          `._  .-~     ~-.   _,'
                           ( )'           '.( )
             `._    _       /               .'
              ( )--' `-.  .'                 ;
         .    .'        '.;                  ()
          `.-.`           '                 .'
----*-----;                                .'
          .`-'.           ,                `.
         '    '.        .';                  ()
              (_)-   .-'  `.                 ;
             ,'   `-'       \               `.
                           (_).           .'(_)
                          .'   '-._   _.-'    `.
                                 .'   `.
                                 '     ; 
                                  `-,-'
                                   /`\
                                 /`
```
