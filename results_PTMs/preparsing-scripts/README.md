# Preparsing scripts

## Coverage

Global coverage (% per histone):

```bash
mkdir -p data-coverage
# also for paired searches?
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/*me2*.pdResult ; do o=$(basename $i) ; echo ${o} ; sqlite3 -header -csv $i < preparsing-scripts/get_coverage_per_protein.sql > data-coverage//${o%%.pdResult}.coverage.csv; done

# get final table
Rscript preparsing-scripts/s20_get_coverage_2021-05-10.R
```

Per-position coverage:

```bash
# get sequence ids
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/*me2*.pdResult ; do o=$(basename $i) ; echo ${o} ; sqlite3 -header -csv $i < preparsing-scripts/get_coverage_per_position-s01.sql > data-coverage//${o%%.pdResult}.perpos_cov_seqids.csv; done

# get coverage per positions
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/*me2*.pdResult ; do o=$(basename $i) ; echo ${o} ; sqlite3 -header -csv $i < preparsing-scripts/get_coverage_per_position-s02.sql > data-coverage//${o%%.pdResult}.perpos_cov_peptides.csv; done
# data tables are named slightly different for Gefoke
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/Gefoke-me2-ac.pdResult ; do o=$(basename $i) ; echo ${o} ; sqlite3 -header -csv $i < preparsing-scripts/get_coverage_per_position-s02-gefoke.sql > data-coverage//${o%%.pdResult}.perpos_cov_peptides.csv; done

# integrate
Rscript preparsing-scripts/s20_get_coverage_per_pos_2021-05-11.R
```

## Paired searches

Data from paired searches of acetylation with mono/di/trimethylation (three files per species). These files are parsed via `sqlite3`.

1. Get CSV tables with modifications from `.pdResult` SQLite databases:

```bash
# general
i=file.pdResult
sqlite3 -header -csv $i < get_PTMRS.sql > ${i%%.pdResult}.csv

# for most species:
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/*.pdResult ; do o=$(basename $i) ; echo -e "\n>>>> $o <<<<" ; sqlite3 -header -csv $i < scripts/get_PTMRS.sql > ../data-paired-searches/${o%%.pdResult}.raw.csv ; done

# for Gefionella specifically:
for i in /users/asebe/xgraubove/histonome-ops/proteome-discoverer-data/data-paired-modifications-search/Gefoke-me*.pdResult ; do o=$(basename $i) ; echo -e "\n>>>> $o <<<<" ; sqlite3 -header -csv $i < scripts/get_PTMRS-fromPsms.sql | sort -u > ../data-paired-searches/${o%%.pdResult}.raw.csv ; done
```

2. Create clean peptide tables:

```bash
Rscript s10_parse_paired_searches_2021-04-26.R

# concatenate them?
for h in H2A macroH2A H2AZ H3 H4 H2B ; do for i in ../data-paired-searches/*-${h}.csv ;  do awk 'NR>1' $i ;  done | sort  -k1,1 -k5,5 -k8,8n -k12,12n > ../data-paired-searches/all_${h}.tsv ; done
```
