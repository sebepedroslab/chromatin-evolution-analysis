# Chromatin evolution

This repository contains all the necessary code to reproduce the comparative proteomics & comparative genomics analyses of our manuscript **A phylogenetic and proteomic reconstruction of eukaryotic chromatin evolution** (Grau-Bové et al., [Nature Ecology and Evolution 2022](https://www.nature.com/articles/s41559-022-01771-6), free PDF available [here](https://www.sebepedroslab.org/_files/ugd/d1ace5_23ff83e58266497f8de5d657ff376708.pdf)).

The data and code are structured as follows:

1. **Classification of canonical and variant histones**: see `results_histones_phylo` folder. Classification of histone variants using phylogenetic and pairwise similarity approaches, and characterisation of N-terminal tails in archaea.

![imatge](https://user-images.githubusercontent.com/11460546/140971995-cf85de81-17ed-4146-a1cf-27f7505430a1.png)

2. **Evolutionary analysis of the chromatin enzymatic machinery and readers**: see `results_toolkit_phylo` folder. Identification of orthologous groups in chromatin modification-related enzyme families from selected eukaryotic species, analyses of their prokaryotic homology, fusions with transposon-associated domains, and more. This folder also includes detailed instructions on how to get all the genomic data used in this study.

![imatge](https://user-images.githubusercontent.com/11460546/140972108-709e020a-a61c-4e3d-a0bc-0fc0c05b6b0f.png)

3. **hPTM conservation**: see `results_PTMs` folder. Includes consensus alignments and our database of histone post-translational modifications with homologous positions in each canonical histone. Here we analyse proteomics data from *Proteome Discoverer* to identify homologous PTMs.

![imatge](https://user-images.githubusercontent.com/11460546/140971814-929ed019-1461-4fef-8748-be4ac8755250.png)

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
