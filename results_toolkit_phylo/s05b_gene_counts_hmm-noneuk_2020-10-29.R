library(data.table)
library(pheatmap)
library(stringr)


# input files
taxne_fn = "../data/taxonomy_ranked_wg.tsv.gz"                 # taxonomic info for non-euks
gen_list_fn = "../data/gene_families_hmm.csv"               # list of gene familes
search_fn = "gene_counts/"
ne_taxon_fo = "../data/"
outp_fn = "results_counts/"

graphics.off()

# heatmap colors
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))
# col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","dodgerblue4"))

# load list hmms
gen_list = read.table(gen_list_fn, header = F, stringsAsFactors = T)
colnames(gen_list) = c("Class","Type","Family","Domains","search_thr","inflation","min_phylo_size")
gen_list_gap_ixs = c(1,1+which(diff(as.numeric(gen_list$Type))!=0)) - 1

#### non-euks HMM counts ####

# create counts tables: for each supergroup (Archaea, Bacteria, Viruses)

# load general taxonomy (species, genus, etc)
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=F, stringsAsFactors=F)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")

noneuk_report_list = list(
  # commented-out taxa have 0 hits
  list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="phylum"),
  list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="family"),
  list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="genus"),
  list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Euryarchaeota", tti="phylum", tfi="genus"),
  list(tax="Archaea", tbi="arc", tai="Euryarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Methanosarcina", tti="genus", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Methanospirillum", tti="genus", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Aenigmarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Altiarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Bathyarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Diapherotrites", tti="phylum", tfi="taxon"),
  # list(tax="Archaea", tbi="arc", tai="Candidatus Geoarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Geothermarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Heimdallarchaeota", tti="phylum", tfi="taxon"),
  # list(tax="Archaea", tbi="arc", tai="Candidatus Helarchaeota", tti="phylum", tfi="taxon"),
  # list(tax="Archaea", tbi="arc", tai="Candidatus Huberarchaea", tti="phylum", tfi="taxon"),
  # list(tax="Archaea", tbi="arc", tai="Candidatus Hydrothermarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Korarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Lokiarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Marsarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Micrarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Nezhaarchaeota", tti="phylum", tfi="taxon"),
  # list(tax="Archaea", tbi="arc", tai="Candidatus Odinarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Parvarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Thorarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Verstraetearchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Candidatus Woesearchaeota", tti="phylum", tfi="taxon"),
  list(tax="Archaea", tbi="arc", tai="Crenarchaeota", tti="phylum", tfi="genus"),
  list(tax="Archaea", tbi="arc", tai="Thaumarchaeota", tti="phylum", tfi="genus"),
  # list(tax="Archaea", tbi="arc", tai="Nanoarchaeota", tti="phylum", tfi="genus"),
  # list(tax="Archaea", tbi="arc", tai="Nanoarchaeota", tti="phylum", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Viruses", tti="superkingdom", tfi="family"),
  list(tax="Viruses", tbi="vir", tai="Nucleocytoviricota", tti="kingdom", tfi="family"),
  list(tax="Viruses", tbi="vir", tai="Caudovirales", tti="order", tfi="family"),
  list(tax="Viruses", tbi="vir", tai="Mimiviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Mimiviridae", tti="family", tfi="genus"),
  list(tax="Viruses", tbi="vir", tai="Marseilleviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Pithoviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Phycodnaviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Iridoviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Myoviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Siphoviridae", tti="family", tfi="taxon"),
  list(tax="Viruses", tbi="vir", tai="Podoviridae", tti="family", tfi="taxon"),
  list(tax="Bacteria", tbi="bac", tai="Bacteria", tti="superkingdom", tfi="phylum"),
  list(tax="Bacteria", tbi="bac", tai="Acidobacteria", tti="phylum", tfi="genus"),
  list(tax="Bacteria", tbi="bac", tai="Bacteroidetes", tti="phylum", tfi="genus"),
  list(tax="Bacteria", tbi="bac", tai="Actinobacteria", tti="phylum", tfi="genus"),
  list(tax="Bacteria", tbi="bac", tai="Firmicutes", tti="phylum", tfi="genus"),
  list(tax="Bacteria", tbi="bac", tai="Spirochaetes", tti="phylum", tfi="genus"),
  list(tax="Bacteria", tbi="bac", tai="Proteobacteria", tti="phylum", tfi="genus")
)

# noneuk_report_list = list(
#   list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="phylum")
# )

# loop through lists
for (rep in noneuk_report_list) {
  
  # define subset to plot
  tax=rep$tax # which file to load? (archaea, virus, bacteria)
  tbi=rep$tbi # brief name for the loaded dataset (arc, vir, bac)
  tai=rep$tai # subset to this taxonomic range
  tti=rep$tti # the subset taxonomic range belongs to this category (kingdom...)
  tfi=rep$tfi # focus on reporting at this lower taxonomic level
  
  print(sprintf("gene counts %s (%s) at %s level", tai, tti, tfi))
  
  # load sequence taxonomy (species of origin of each seq)
  gentax_fn = sprintf("%s/seq_%s.taxa.csv.gz", ne_taxon_fo, tax)
  gentax = fread(cmd=sprintf("zcat %s", gentax_fn), data.table = F, header = F)
  colnames(gentax) = c("gene", "taxon")
  
  # restrict taxonomy to group of interest
  taxni = taxne[!is.na(taxne[,tti]),]
  taxni = taxni[taxni[,tti] == tai,]
  # restrict taxonomy to taxa that actually appear on the sequence set
  taxni = taxni [ taxni$taxon %in% gentax$taxon, ]
  # taxni = droplevels(taxni)
  taxni_counts = table(taxni[,tfi])
  taxni_counts = c(sum(taxni_counts),taxni_counts)
  names(taxni_counts)[1] = "Sum"
  print(sprintf("gene counts across %i %s %s(s)", length(unique(taxni[,tfi])), tai, tfi))
  
  # load HMM hits
  gen = fread(input = sprintf("gene_counts/%s_genecounts.csv", tbi), data.table = F, header = F, col.names = c("gene","gene_family"))
  
  # add species info
  gen = merge(gen, gentax, by.x = "gene", by.y = "gene", all.x = T)
  # one row per species
  gen_u = within(gen, rm("gene"))
  gen_u = unique(gen_u)

  # add taxonomic info
  gen_wtax = merge(gen_u, taxni, by.x = "taxon", by.y = "taxon", all.x = T, all.y = F)
  
  # taxonomic groupings and gene families as factors
  if (tax=="Archaea" & tbi == "arc" &  tfi == "phylum" & tti == "superkingdom") {
    levels_taxon_factor = c(
      "Candidatus Lokiarchaeota","Candidatus Thorarchaeota","Candidatus Helarchaeota","Candidatus Odinarchaeota",
      "Candidatus Heimdallarchaeota","Crenarchaeota","Candidatus Marsarchaeota","Candidatus Verstraetearchaeota",
      "Candidatus Nezhaarchaeota","Candidatus Bathyarchaeota","Thaumarchaeota","Aigarchaeota",
      "Candidatus Geothermarchaeota","Candidatus Korarchaeota","Euryarchaeota","Candidatus Diapherotrites",
      "Candidatus Micrarchaeota","Candidatus Altiarchaeota","Candidatus Hydrothermarchaeota","Candidatus Aenigmarchaeota",
      "Candidatus Nanohaloarchaeota",
      "Candidatus Huberarchaea","Candidatus Parvarchaeota","Nanoarchaeota","Candidatus Woesearchaeota")
  } else {
    levels_taxon_factor = levels(as.factor(taxni[,tfi]))
  }
  gen_wtax[,"taxon_factor"] = factor(gen_wtax[,tfi], levels = levels_taxon_factor)
  gen_wtax[,"genefam_factor"] = factor(gen_wtax$gene_family, levels = as.vector(gen_list$Family))
  
  # crosstabulation
  gen_crosstab = xtabs(formula = ~ genefam_factor + taxon_factor, data = gen_wtax, drop.unused.levels	= F)
  if (tax=="Archaea" & tbi == "arc" &  tfi == "phylum" & tti == "superkingdom") {
    gen_crosstab_tax_present = apply(gen_crosstab, 2, sum) >= 0
  } else {
    gen_crosstab_tax_present = apply(gen_crosstab, 2, sum) > 0
  }
  gen_crosstab = gen_crosstab[,gen_crosstab_tax_present]
  gen_crosstab = cbind(rowSums(gen_crosstab, na.rm = T), gen_crosstab)
  colnames(gen_crosstab)[1] = "Sum"
  gen_crosstab_row = rownames(gen_crosstab)
  gen_crosstab_col = colnames(gen_crosstab)
  gen_crosstab = as.data.frame(t(matrix(gen_crosstab, nrow = nrow(gen_crosstab))))
  colnames(gen_crosstab) = gen_crosstab_row
  rownames(gen_crosstab) = gen_crosstab_col

  # plot absolute count (number of species with gene presence within group)
  # define table splits (for very large tables)  
  mat_pres_cla_split = split(gen_crosstab, rep(1:ceiling(nrow(gen_crosstab)/100), each=100, length.out=nrow(gen_crosstab)))
  pdf(file=sprintf("%s/counts_noneuks_genes_%s_%s.pdf", outp_fn, gsub(" ","_", x=tai), tfi),height=10,width=8)
  non=0
  for (mat_pres_cla_i in mat_pres_cla_split) {
    
    non=non+1
    pheatmap(mat_pres_cla_i, 
             color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
             gaps_col = gen_list_gap_ixs,
             gaps_row = 1,
             cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_format = "%i",
             main=sprintf("Gene presence per %s in %s (%s) (%i/%i)", tfi, tai, tti,non, length(mat_pres_cla_split)))
    
  }
  
  # now plot fraction of species within group  
  # define table splits (for very large tables)  
  gen_crosstab_frac = sweep(gen_crosstab, MARGIN=1, taxni_counts[rownames(gen_crosstab)], "/")
  mat_pres_cla_split = split(gen_crosstab_frac, rep(1:ceiling(nrow(gen_crosstab)/100), each=100, length.out=nrow(gen_crosstab)))
  non=0
  for (mat_pres_cla_i in mat_pres_cla_split) {
    
    non=non+1
    pheatmap(mat_pres_cla_i, 
             color = col_blue(20), breaks = c(0,0.0001,seq(0.05,1,length.out = 19)), 
             gaps_col = gen_list_gap_ixs,
             gaps_row = 1,
             cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = F, number_format = "%i",
             main=sprintf("Gene presence per %s in %s (%s) (%i/%i) (fraction of species)", tfi, tai, tti,non, length(mat_pres_cla_split)))
    
  }
  
  # plot num sps per group
  par(mar=c(40, 4.1, 4.1, 2.1))
  barplot(taxni_counts[rownames(gen_crosstab)], 
          las=2, border = NA, cex.names = 0.5,
          names.arg = paste(names(taxni_counts[rownames(gen_crosstab)]),"n =", taxni_counts[rownames(gen_crosstab)]),
          main="num species per group")
  
  dev.off()
  
  # # plot counts per group
  # pheatmap(gen_crosstab, 
  #          color = col_blue(10), breaks = seq(0, 10, length.out = 11)-0.01, 
  #          gaps_col = gen_list_gap_ixs,
  #          cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5,
  #          border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_format = "%i",
  #          main=sprintf("Gene presence per %s in %s (%s)", tfi, tai, tti))
  
  # save csvs
  gen_crosstab = cbind(gen_crosstab, data.frame(total = taxni_counts[rownames(gen_crosstab)]))
  gen_crosstab = within(gen_crosstab, rm("total.Var1"))
  write.table(
    gen_crosstab, file=sprintf("%s/counts_noneuks_genes_%s_%s.csv", outp_fn, gsub(" ","_", x=tai), tfi), 
    sep="\t",
    quote = F)
  
}



print("Done!")