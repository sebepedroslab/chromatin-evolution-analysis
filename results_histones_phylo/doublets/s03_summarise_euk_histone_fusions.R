# libraries
library(stringr)
library(GenomicRanges)
library(pheatmap)
library(ape)


# species newick
phyl_fn = "../../results_toolkit_phylo/data/species_tree.newick"


# load dimers table
hdi = read.table(file="histone_dimers.csv", sep="\t", stringsAsFactors = FALSE, header = TRUE)

# fetch euks
hdie = hdi[is.na(hdi$taxon),c("gene_i","classification_i","classification_j","start_i","end_i","start_j","end_j")]
hdie$fusion = paste(hdie$classification_i, hdie$classification_j)
hdie$species = stringr::str_split(hdie$gene_i, pattern = "_", simplify = TRUE)[,1]
colnames(hdie)[colnames(hdie) == "gene_i"] = "gene"

# for each eukaryotic species, check fusion completeness diagnostics
sps_list = unique(hdie$species)
sps_list = c("Azfi","Ddis","Dpul","Drer","Hetalb","Klenit","Mgut","Nvec","Perkma","Ppat",
             "Sphfal","Symmic","Ttra","Apla","Spur","Adig","Exapal","Skow","Ctel","Spis")
gen_fusions = data.frame()

# order species (from phylogeny)
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label

sps_list = sps_list[match(sps_order, sps_list)]
sps_list = sps_list[!is.na(sps_list)]

for (sps in sps_list) {
  
  print(sps)
  
  #### Define per-species candidates ####
  
  # create dataframe to store per-gene diagnostics
  # spi_fusions = data.frame(gene = hdie[hdie$species == sps,"gene_i"], architecture = hdie[hdie$species == sps,"fusion"])
  spi_fusions = hdie[hdie$species == sps,]
  rownames(spi_fusions) = spi_fusions$gene
  
  # read pfam domain coordinates
  spi_pfam = read.table(sprintf("data/fusions.%s.Pfamscan.seqs.csv", sps), stringsAsFactors = FALSE,
                        header = FALSE, col.names = c("gene","dom_start","dom_end","pfam_id","pfam_dom","pfam_seq"))
  
  
  #### Assembly contiguity ####
  
  # read assembly contiguity
  # transcripts that are contiguous in the genome assembly (no NNN stretches within gene model)
  if (file.info(sprintf("data/fusions.%s.nonNstretch.txt", sps))$size > 0) {
    spi_acon = read.table(sprintf("data/fusions.%s.nonNstretch.txt", sps), stringsAsFactors = FALSE,
                          header = FALSE)[,1]
    spi_fusions$is_contiguous_assembly = spi_fusions$gene %in% spi_acon
  } else {
    # if there is no assembly contiguity data, record as FALSE
    spi_fusions$is_contiguous_assembly = FALSE
  }
  
  
  #### Expression ####
  
  # read expression
  spi_expr = read.table(sprintf(
    "data/fusions.%s.salmon_quant.csv", sps), 
    stringsAsFactors = FALSE, header = TRUE)
  if (nrow(spi_expr) > 1) { 
    spi_expr_s = spi_expr[c("Gene","TPM","Sample")]
    spi_expr_m = reshape(spi_expr_s, idvar = "Gene", timevar = "Sample", direction = "wide")
    rownames(spi_expr_m) = spi_expr_m$Gene
    spi_expr_m = within(spi_expr_m, rm("Gene"))
    colnames(spi_expr_m) = stringr::str_split(colnames(spi_expr_m), pattern = "\\.", simplify = TRUE)[,2]
    # number of samples in which each gene is expressed
    spi_expr_count = data.frame(gene = names(rowSums(spi_expr_m > 0)), num_expressed=rowSums(spi_expr_m > 0))
    spi_fusions = merge(spi_fusions, spi_expr_count, all.x = TRUE, all.y = FALSE)
    spi_fusions$num_expressed [is.na(spi_fusions$num_expressed )] = 0
    spi_fusions$is_expressed = spi_fusions$num_expressed > 0
    # read expression contiguity
    spi_exco = read.table(sprintf("data/fusions.%s.salmon_cover.csv", sps), header = TRUE, stringsAsFactors = FALSE)
    spi_exco = merge(spi_exco, spi_expr, by.x="Transcript", by.y="Transcript")
    spi_exco_gr = GenomicRanges::makeGRangesFromDataFrame(spi_exco, start.field = "start_pep", end.field = "end_pep", seqnames.field = "Gene")
    spi_full_coverage = c() # empty vector (gets filled below)
    
    
    #### Contiguous expression coverage ####
    
    # pfam coordinates per gene to identify core and TE region
    # - core region: first core domain found in the gene
    # - TE region: last TE domain found in the gene
    # - result: the gap will be one of the most extensive in that gene that bridges the fusion
    # define GRanges of core region
    spi_pcor_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_fusions, start.field = "start_i",end.field = "end_i", seqnames.field = "gene")
    # define GRanges of TE region
    spi_ptes_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_fusions, start.field = "start_j",end.field = "end_j", seqnames.field = "gene")
    # define GRanges gaps between core and TE
    spi_pgap_gr = IRanges::pgap(spi_pcor_gr, spi_ptes_gr)
    
    # overlap between domains and coverage deserts
    spi_ovp_pgap_exco = IRanges::findOverlaps(query = spi_pgap_gr, spi_exco_gr, type = "within")
    spi_full_coverage = c(spi_full_coverage,as.vector(spi_pgap_gr@seqnames)[unique(spi_ovp_pgap_exco@from)])
    
    # report full coverage
    spi_fusions$is_contiguous_coverage = spi_fusions$gene %in% spi_full_coverage
    
  } else {
    
    # if there is no expression data, record as FALSE
    spi_fusions$num_expressed = NA
    spi_fusions$is_expressed = FALSE
    spi_fusions$is_contiguous_coverage = FALSE
    
  }
  
  #### Store data ####
  gen_fusions = rbind(gen_fusions, spi_fusions)
  
}


#### Summary #####

# species order
gen_fusions$species = factor(gen_fusions$species, levels = sps_list)

# summarise expression status
gen_fusions$evidence_status = "no evidence"
gen_fusions[gen_fusions$is_expressed  & gen_fusions$is_contiguous_coverage  & gen_fusions$is_contiguous_assembly,"evidence_status"]  = "contiguous assembly, complete expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_coverage & gen_fusions$is_contiguous_assembly,"evidence_status"]  = "contiguous assembly, partial expression"
gen_fusions[!gen_fusions$is_expressed & !gen_fusions$is_contiguous_coverage & gen_fusions$is_contiguous_assembly,"evidence_status"]  = "contiguous assembly, no expression"
gen_fusions[!gen_fusions$is_expressed & gen_fusions$is_contiguous_coverage  & gen_fusions$is_contiguous_assembly,"evidence_status"]  = "contiguous assembly, no expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_assembly & gen_fusions$is_contiguous_assembly,"evidence_status"]  = "broken assembly, complete expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_assembly & !gen_fusions$is_contiguous_assembly,"evidence_status"] = "broken assembly, partial expression"
gen_fusions[!gen_fusions$is_expressed & !gen_fusions$is_contiguous_assembly,"evidence_status"] = "broken assembly, no expression"
# factors
gen_fusions$evidence_status = factor(
  gen_fusions$evidence_status,
  levels = c("contiguous assembly, complete expression",
             "contiguous assembly, partial expression",
             "contiguous assembly, no expression",
             "broken assembly, complete expression",
             "broken assembly, partial expression",
             "broken assembly, no expression",
             "no evidence"))

# save table
pdf("gene_fusion_evidence.pdf", height = 5, width = 5)
par(mar = c(5.1, 20, 4.1, 2.1))
barplot(
  table(gen_fusions$evidence_status), horiz = TRUE, las=1,
  names.arg = paste(names(table(gen_fusions$evidence_status)),"n =", table(gen_fusions$evidence_status)),
  main = sprintf("histone fusions\nn=%i", nrow(gen_fusions), cex.names = 0.6, cex.axis = 0.6)
)

barplot(
  table(gen_fusions$fusion), horiz = TRUE, las=1,
  names.arg = paste(names(table(gen_fusions$fusion)),"n =", table(gen_fusions$fusion)),
  main = sprintf("histone fusions\nn=%i", nrow(gen_fusions), cex.names = 0.3, cex.axis = 0.6)
)
dev.off()




## heatmaps counts
col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","midnightblue"))

pdf("heatmap_histone_fusion_evidence.pdf", height = 6, width = 5)


#### evidence + fusion type ####
tab = xtabs(formula = ~ evidence_status + fusion, data = gen_fusions, drop.unused.levels	= TRUE)
tam = matrix(tab, nrow = dim(tab)[1])
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
pheatmap(tam [ , order(colSums(tam), decreasing = TRUE) ],  
         color = col_blue(10), breaks = seq(0,10,length.out = 11) - 0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize = 5, 
         main = sprintf("summary fusions"),
         border_color = "white", display_numbers = TRUE, number_format = "%i")


#### evidence + fusion type + species ####
gen_fusions_i = gen_fusions

# order rows according to combinations of species and core domains (proxy to gene family)
gen_fusions_i$species = droplevels(gen_fusions_i$species)
# keep species
keep_species = unique(as.character(gen_fusions_i$species))

# col separation
colsep = c()

# table TE classification
tab = xtabs(formula = ~ species + fusion, data = gen_fusions_i, drop.unused.levels	= FALSE)
tam = matrix(tab, nrow = dim(tab)[1])
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
tat = tam [ , order(colSums(tam), decreasing = TRUE) ]
colsep = c(colsep, ncol(tat))

# table evidence status
tab = xtabs(formula = ~ species + evidence_status, data = gen_fusions_i, drop.unused.levels	= FALSE)
tam = matrix(tab, nrow = dim(tab)[1])
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
tat = cbind(tat, tam)
colsep = c(colsep, ncol(tat))

# plot
tas = matrix(tat[keep_species,], nrow = length(keep_species))
rownames(tas) = keep_species
colnames(tas) = colnames(tat)
pheatmap(tas,  
         gaps_col = colsep,
         color = col_blue(10), breaks = seq(0,10,length.out = 11) - 0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize = 5, 
         main = sprintf("all fusions"),
         border_color = "white", display_numbers = TRUE, number_format = "%i")

#### evidence + fusion type + species, per type ####
for (cli in unique(gen_fusions$fusion)) {
  
  # subset
  gen_fusions_i = gen_fusions[gen_fusions$fusion == cli,]
  
  # order rows according to combinations of species and core domains (proxy to gene family)
  gen_fusions_i$species = droplevels(gen_fusions_i$species)
  # keep species
  keep_species = unique(as.character(gen_fusions_i$species))
  
  
  if (nrow(gen_fusions_i) > 0) {
    # col separation
    colsep = c()
    
    # table TE classification
    tab = xtabs(formula = ~ species + fusion, data = gen_fusions_i, drop.unused.levels	= FALSE)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = tam
    colsep = c(colsep, ncol(tat))
    
    # table evidence status
    tab = xtabs(formula = ~ species + evidence_status, data = gen_fusions_i, drop.unused.levels	= FALSE)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = cbind(tat, tam)
    colsep = c(colsep, ncol(tat))
    
    # plot
    tas = matrix(tat[keep_species,], nrow = length(keep_species))
    rownames(tas) = keep_species
    colnames(tas) = colnames(tat)
    pheatmap(tas,  
             gaps_col = colsep,
             color = col_blue(10), breaks = seq(0,10,length.out = 11) - 0.01, 
             cellwidth = 4, cellheight = 4, na_col = "grey", 
             cluster_rows = FALSE, cluster_cols = FALSE,
             fontsize = 5, 
             main = sprintf("%s fusions", cli),
             border_color = "white", display_numbers = TRUE, number_format = "%i")
  }
}


dev.off()







# write output
write.table(
  gen_fusions, 
  file = "histone_dimers_evidence.csv" ,
  quote = FALSE, sep = "\t", row.names = FALSE)

# summary of evidence per fusion type and species:
options(width = 1000) 
capture.output(
  xtabs(formula = ~ evidence_status + species + fusion, data = gen_fusions, drop.unused.levels	= TRUE),
  file = "histone_dimers_evidence_crosstab.txt")

print("Done!")



#### fasta output ####

fas = seqinr::read.fasta("fusions_all.fasta")

keep_seq = gen_fusions[gen_fusions$evidence_status == "contiguous assembly, complete expression", "gene"]
keep_fus = gen_fusions[gen_fusions$evidence_status == "contiguous assembly, complete expression", "fusion"]
keep_fus = gsub(" ","_", keep_fus)
keep_fus = gsub("other_","", keep_fus)

fai = fas[keep_seq]

for (i in 1:length(fai)) {
  fi = fai[i]
  fu = keep_fus[i]
  seqinr::write.fasta(sequences = fi[[1]], names = names(fi), file.out = sprintf("valid-fusions/%s_%s.fasta", names(fi), fu))
}
