# Alignments

1. MSA from `mafft`:

```bash
mafft --globalpair --thread 1 --reorder --maxiterate 1000 input.fasta > output.fasta
```

2. Inspect alignment & select consensus sequence (keep close to human orthologs).

3. Remove initial `M` residue from consensus alignment (hPTM coordinates are given from first non-`M` residue).

4. Realign.

## Archaea

* `Methanospirillum_stamsii_WP_109941938.1|Histone`: seems closer to H3, few homologous positions (maybe H3K4?).

* `Methanosarcina_spelaei_PAV14247.1|Histone`: seems closer to H4, few homologous positions.

* Other histones
