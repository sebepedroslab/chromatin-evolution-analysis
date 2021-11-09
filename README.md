# Chromatin evolution

This repository contains all the necessary code to reproduce the comparative proteomics & comparative genomics analyses of our manuscript [**Comparative proteogenomics deciphers the origin and evolution of eukaryotic chromatin**](https://github.com/sebepedroslab/chromatin-evolution-analysis) (Grau-Bové et al. 2021).

The data and code are structured as follows:

1. **Classification of canonical and variant histones**: see `results_histones_phylo` folder. Classification of histone variants using phylogenetic and pairwise similarity approaches, and characterisation of N-terminal tails in archaea.

2. **Evolutionary analysis of the chromatin enzymatic machinery and readers**: see `results_toolkit_phylo` folder. Identification of orthologous groups in chromatin modification-related enzyme families from selected eukaryotic species, analyses of their prokaryotic homology, fusions with transposon-associated domains, and more. This folder also includes detailed instructions on how to get all the genomic data used in this study.

3. **hPTM conservation**: see `results_PTMs` folder. Includes consensus alignments and our database of histone post-translational modifications with homologous positions in each canonical histone. Here we analyse proteomics data from *Proteome Discoverer* to identify homologous PTMs.

```python
 .              +   .                .   . .     .  .
                   .                    .       .     *
  .       *                        . . . .  .   .  + .
          Hic sunt dracones            .   .  +  . . .
.                 |             .  .   .    .    . .
                  |           .     .     . +.    +  .
                 \|/            .       .   . .
        . .       V          .    * . . .  .  +   .
           +      .           .   .      +
                            .       . +  .+. .
  .                      .     . + .  . .     .      .
           .      .    .     . .   . . .        ! /
      *             .    . .  +    .  .       - O -
          .     .    .  +   . .  *  .       . / |
               . + .  .  .  .. +  .
.      .  .  .  *   .  *  . +..  .            *
 .      .   . .   .   .   . .  +   .    .            +
```
