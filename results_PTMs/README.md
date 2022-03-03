# Data folders

Data & data folders:

* **`consensus_modifications_perseq.xlsx`**: summary of homologous PTM positions in eukaryotes (parsed from ProteomeDiscoverer SQL datasets).
* **`species_list_2021-02-10.txt`**: list of eukaryotic species covered in this analysis.
* `data-coverage/`: per-position and per-sequence coverage data for the paired searches.
* `data-archaea/`: alignments of archaeal histones for four species (see readme inside).

Scripts:

```bash
# create heatmaps and summary tables with per-species and per-position counts of PTM evidence
Rscript s01_plot_hPTM_tables_2020-10-08.R
Rscript s02_plot_ubi_tables_2021-03-08.R

# create summaries of per-sequence and per-position proteomics coverage data
Rscript s20_get_coverage_2021-05-10.R
Rscript s21_get_coverage_per_pos_2021-05-11.R

# evaluate reproducibility of histone peptide-spectrum matches with hPTMs:
Rscript s22_get_replicability_per_peptide_2022-02-16.R
```
