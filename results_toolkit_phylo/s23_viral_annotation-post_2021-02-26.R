# load libraries
library("viridis")
library("scales")
library("stringr")
library("pheatmap")
library("plotrix")

#### Define input ####

virt = read.table("results_viruses/summary_viral_homologs_annotation.csv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# domains with viral hits
domlist=c("SET","BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","HIF-1","Histone","YEATS",
          "LinkerHistone","PHD","PTIP","SIR2","SNF2_N","TUDOR",
          "zf-CCHH","zf-CXXC","Hist_deacetyl","Acetyltransf_1")

# tables
ort_fn = "orthogroups_euk.csv"
cla_fn = "../data/gene_families_hmm.csv"
arq_fn = "gene_sequences/noneuk_vir.pfamscan_archs.csv"

# read data
# orthology, gene fam classification, architectures
ort = read.table(ort_fn, sep="\t", header = TRUE, stringsAsFactors = FALSE)
cla = read.table(
  cla_fn, header = FALSE, sep = "\t", 
  col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
  stringsAsFactors = FALSE)
arq = read.table(arq_fn, header = FALSE, col.names = c("gene","architecture"), sep="\t")

# simplify architectures (collapse **consecutive** repeats)
simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
arq$simarq = simple_arqs_vect

# ordered list of domains
domlist = domlist[match(cla$gene_fam, domlist)]
domlist = domlist[!is.na(domlist)]

# get gene class list
clalist = unique(cla$gene_type)

# get families of relevant groups, to create appropriate factor
# first, large DNA viruses
virt_giantdna_fams = unique(virt[virt$kingdom == "Nucleocytoviricota",]$family)
virt_giantdna_fams = sort(virt_giantdna_fams[!is.na(virt_giantdna_fams)])
# second, Caudovirales
virt_caudovir_fams = unique(virt[virt$order == "Caudovirales",]$family)
virt_caudovir_fams = sort(virt_caudovir_fams[!is.na(virt_caudovir_fams)])
# finally, the rest
virt_othertax_fams = unique(virt$family)
virt_othertax_fams = sort(virt_othertax_fams[!is.na(virt_othertax_fams)])
virt_othertax_fams = virt_othertax_fams[ !virt_othertax_fams %in% virt_giantdna_fams ]
virt_othertax_fams = virt_othertax_fams[ !virt_othertax_fams %in% virt_caudovir_fams ]
# list of taxa groups
tax_list = c(virt_giantdna_fams, virt_caudovir_fams , virt_othertax_fams, "Unknown")
# now create factor
virt$tax_factor = factor(
  virt$family,
  levels = tax_list)
virt$tax_factor[is.na(virt$family)] = "Unknown"
# colors
tax_list_col = rainbow(length(tax_list), end = 0.9, v=0.9)

# supra taxon
virt$tax_factor_supra = "Other"
virt$tax_factor_supra [ virt$tax_factor %in% virt_giantdna_fams ] = "Nucleocytoviricota"
virt$tax_factor_supra [ virt$tax_factor %in% virt_caudovir_fams ] = "Caudovirales"
virt$tax_factor_supra = factor(virt$tax_factor_supra, levels = c("Nucleocytoviricota", "Caudovirales","Other"))


# factor of gene families
virt$core_domain_factor = factor(virt$core_domain, levels = domlist)
# color
dom_list_col = rainbow(length(domlist), end = 0.9, v=0.9)

# factor for closest group
virt$closest_group_factor = factor(virt$closest_group_filt, levels = c("euk","arc","bac","unk"))
virt$closest_group_factor[is.na(virt$closest_group)] = "unk"
cgo_list_col = c("chartreuse3","slategray4","lightgray","darkolivegreen2")



# barplots
# summary counts
pdf("results_viruses/summary_viral_homologs_barplots.pdf", height = 8, width = 6)
par(mar = c(5.1, 16, 4.1, 2.1))

# taxa per dom
xtab = table(virt$tax_factor,virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
barplot(xtab, col = tax_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32))
legend("bottomright", fill=tax_list_col, legend = tax_list, cex=0.6, bty="n")
title(sprintf("viral taxa across domains\nn=%i",sum(xtab)))

# taxa per dom
xtab = table(virt$tax_factor,virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
xtab = t(xtab)
xtab_frac = xtab / rowSums(xtab)
xtab_frac = t(xtab_frac)
barplot(xtab_frac, col = tax_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32), xlim = c(0,1))
# legend("bottomright", fill=tax_list_col, legend = tax_list, cex=0.6, bty="n")
title(sprintf("viral taxa across domains, fraction\nn=%i",sum(xtab)))

# closest annotation per dom
xtab = table(virt$closest_group_factor, virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
barplot(xtab, col = cgo_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32))
legend("bottomright", fill=cgo_list_col, legend = c("euk","arc","bac","unk"), cex=0.6, bty="n")
title(sprintf("closest annotation across domains\nn=%i",sum(xtab)))

# domain per taxa
xtab = table(virt$core_domain_factor, virt$tax_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
barplot(xtab, col = dom_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32))
legend("bottomright", fill=dom_list_col, legend = domlist, cex=0.6, bty="n")
title(sprintf("domains across taxa\nn=%i",sum(xtab)))

# closest annotation per taxa
xtab = table(virt$closest_group_factor, virt$tax_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
barplot(xtab, col = cgo_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32))
legend("bottomright", fill=cgo_list_col, legend = c("euk","arc","bac","unk"), cex=0.6, bty="n")
title(sprintf("closest annotation across taxa\nn=%i",sum(xtab)))


# closest annotation per core domain
xtab = table(virt$closest_group_factor, virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
barplot(xtab, col = cgo_list_col, horiz = TRUE, las=1, names.arg = xtab_name, ylim = c(0,32))
legend("bottomright", fill=cgo_list_col, legend = c("euk","arc","bac","unk"), cex=0.6, bty="n")
title(sprintf("closest annotation across core domains\nn=%i",sum(xtab)))

# closest annotation per core domain, fraction
xtab = table(virt$closest_group_factor, virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
xtab = t(xtab)
xtab_frac = xtab / rowSums(xtab)
xtab_frac = t(xtab_frac)
barplot(xtab_frac, col = cgo_list_col, horiz = TRUE, las=1, names.arg = xtab_name, xlim = c(0,1), ylim = c(0,32))
legend("bottomright", fill=cgo_list_col, legend = c("euk","arc","bac","unk"), cex=0.6, bty="n")
title(sprintf("closest annotation across core domains\nn=%i",sum(xtab)))

# supra taxa per core domain, fraction
xtab = table(virt$tax_factor_supra, virt$core_domain_factor)
xtab_name = paste(colnames(xtab), "n =", colSums(xtab))
xtab = t(xtab)
xtab_frac = xtab / rowSums(xtab)
xtab_frac = t(xtab_frac)
barplot(xtab_frac, col = c("red","orange","gray"), horiz = TRUE, las=1, names.arg = xtab_name, xlim = c(0,1), ylim = c(0,32))
legend("bottomright", fill=c("red","orange","gray"), legend = c("Giant","Caudovirales","other"), cex=0.6, bty="n")
title(sprintf("closest annotation across supra taxa\nn=%i",sum(xtab)))


dev.off()


# heatmaps per gene family
# cross tabulations of domain annotations per family and taxonomic group
# also add best orthogroup???
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))
tax_gaps = c(length(virt_giantdna_fams),
             length(virt_giantdna_fams) + length(virt_caudovir_fams),
             length(virt_giantdna_fams) + length(virt_caudovir_fams) + length(virt_othertax_fams))

# loop
pdf("results_viruses/summary_viral_homologs_tables.pdf", height = 6, width = 5)
for (dom in domlist) {
  
  # subset  
  virt_i = virt[virt$core_domain == dom,]
  
  # tables of closest homologs per taxa
  xtab1 = table(virt_i$tax_factor, virt_i$closest_group_factor)
  # table of closest orthogroup (side by side)
  xtab2 = table(virt_i$tax_factor, as.character(virt_i$closest_OG_euk_filt))
  colnames(xtab2) = str_trunc(colnames(xtab2), 40)
  xtab2 = xtab2[,order(colSums(xtab2), decreasing = TRUE)]
  # plot side by side
  xtab = cbind(xtab1, xtab2)
  pheatmap(xtab, 
           color = col_blue(10), breaks = seq(0,10,length.out = 11) - 0.01, 
           gaps_row = tax_gaps, gaps_col = ncol(xtab1),
           cellwidth = 5, cellheight = 5, na_col = "white",number_color = "aliceblue", fontsize = 5,
           border_color = "white", cluster_cols=FALSE, cluster_rows=FALSE,display_numbers = TRUE, number_format = "%i",
           main=sprintf("%s annotated OGs, per taxa\nn=%i viral genes", dom,nrow(virt_i)))
  
#   # tables of closest homologs per taxa
#   xtab1 = table(virt_i$tax_factor_supra, virt_i$closest_group_factor)
#   # xtab1 = xtab1/rowSums(xtab1)
#   pheatmap(xtab1, 
#            color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
#            cellwidth = 5, cellheight = 5, na_col = "white",number_color = "aliceblue", fontsize = 5,
#            border_color = "white", cluster_cols=FALSE, cluster_rows=FALSE,display_numbers = TRUE,   number_format = "%i",
#            main=sprintf("%s annotated OGs, per taxa\nn=%i viral genes", dom,nrow(virt_i)))
  
#   # sum fraction of sequence neighbours from each cellular group
#   atab = data.frame(
#     row = aggregate(virt_i$euk_top_frac, by=list(virt_i$tax_factor), function(x) sum(x, na.rm = TRUE)/length(x))[,1],
#     euk = aggregate(virt_i$euk_top_frac, by=list(virt_i$tax_factor), function(x) sum(x, na.rm = TRUE)/length(x))[,2],
#     arc = aggregate(virt_i$arc_top_frac, by=list(virt_i$tax_factor), function(x) sum(x, na.rm = TRUE)/length(x))[,2],
#     bac = aggregate(virt_i$bac_top_frac, by=list(virt_i$tax_factor), function(x) sum(x, na.rm = TRUE)/length(x))[,2]
#   )
#   atab = merge(atab, data.frame(tax=levels(virt_i$tax_factor)), by.x = "row", by.y = "tax", all.x = TRUE, all.y = TRUE)
#   rownames(atab) = atab$row
#   atab = atab[levels(virt_i$tax_factor),]
#   atab = atab[,c("euk","arc","bac")]
#   # table
#   pheatmap(atab,
#            color = col_blue(20), breaks = seq(0,1,length.out = 21), 
#            gaps_row = tax_gaps, 
#            cellwidth = 5, cellheight = 5, na_col = "gray90",number_color = "aliceblue", fontsize = 5,
#            border_color = "white", cluster_cols=FALSE, cluster_rows=FALSE,display_numbers = FALSE, number_format = "%.2f",
#            main=sprintf("%s cum fraction of neighbours, per taxa\nn=%i viral genes", dom,nrow(virt_i)))
  
}
dev.off()




# architectures in viral sequences from each domain
pdf("results_viruses/summary_viral_homologs_tables-arch.pdf", height = 6, width = 5)

for (cli in clalist) {
  
  # subset  
  dom = cla[cla$gene_type == cli, ]$gene_fam
  virt_i = virt[virt$core_domain %in% dom,]
  
  if (any(dom %in% domlist)) {
    
    # cross tab  
    xtab = table(paste(virt_i$core_domain,virt_i$simarq,sep="|"), virt_i$closest_group_factor)
    xtab = xtab[order(rowSums(xtab), decreasing = TRUE),]
    if (is.null(nrow(xtab))) { 
      xtab = t(xtab) 
      rownames(xtab) = unique(virt_i$simarq)
    }
    pheatmap(xtab, 
             color = col_blue(10), breaks = seq(0,10,length.out = 11) - 0.01, 
             cellwidth = 5, cellheight = 5, na_col = "white",number_color = "aliceblue", fontsize = 5,
             labels_row = str_trunc(rownames(xtab), 60),
             border_color = "white", cluster_cols=FALSE, cluster_rows=FALSE,display_numbers = TRUE, number_format = "%i",
             main=sprintf("%s architectures, per closest taxon\nn=%i viral genes", dom,nrow(virt_i)))
    
    # fraction within lineage
    xtab_f = apply(xtab, 1, function(x) x / sum(x, na.rm = TRUE))
    xtab_f [ is.na(xtab_f) ] = 0
    if (is.null(nrow(xtab_f))) { 
      xtab_f = t(xtab_f)
      rownames(xtab_f) = unique(virt_i$simarq)
    }
    pheatmap(t(xtab_f), 
             color = col_blue(20), breaks = seq(0,1,length.out = 21), 
             cellwidth = 5, cellheight = 5, na_col = "white",number_color = "aliceblue", fontsize = 5,
             labels_row = str_trunc(rownames(xtab), 60),
             border_color = "white", cluster_cols=FALSE, cluster_rows=FALSE,display_numbers = TRUE, number_format = "%.2f",
             main=sprintf("%s architectures, per closest taxon\nn=%i viral genes", dom,nrow(virt_i)))
  }
}
dev.off()






# closest homology per architecture
pdf("results_viruses/summary_viral_homologs_tables-arch-closest.pdf", height = 7, width = 12)
layout(mat=matrix(1:9,nrow=3))
par(mar = c(5.1, 20, 4.1, 2.1))
for (cli in clalist) {
  
  # subset  
  dom = cla[cla$gene_type == cli, ]$gene_fam
  virt_i = virt[virt$core_domain %in% dom,]
  
  if (any(dom %in% domlist)) {
    # cross tab  
    xtab = table(paste(virt_i$core_domain,virt_i$simarq,sep=" | "), virt_i$closest_group_factor)
    xtab = xtab[order(rowSums(xtab), decreasing = TRUE),]
    if (is.null(nrow(xtab))) { 
      xtab = t(xtab) 
      rownames(xtab) = unique(virt_i$simarq)
    }
    xtab = head(xtab,16)
    barplot(
      t(xtab), 
      col = c("chartreuse3","slategray4","lightgray","darkolivegreen2"),
      horiz = TRUE,
      names.arg = paste(str_trunc(rownames(xtab),60), " n=",rowSums(xtab),sep = ""  ),
      las=1, cex.names = 0.6, cex.axis = 0.6, border = NA,ylim = c(0,16), xlim = c(0,160))
    legend(
      "topright", 
      fill = c("chartreuse3","slategray4","lightgray","darkolivegreen2"),
      legend = c("euk","arc","bac","unk"), bty="n", cex = 0.3, border = NA)
    title(sprintf("%s\nn=%i",cli, sum(xtab)), cex.main=0.8)
  }
  
}
dev.off()

# taxonomic family per architecture
pdf("results_viruses/summary_viral_homologs_tables-arch-virfam.pdf", height = 7, width = 12)
layout(mat=matrix(1:9,nrow=3))
par(mar = c(5.1, 20, 4.1, 2.1))
for (cli in clalist) {
  
  # subset  
  dom = cla[cla$gene_type == cli, ]$gene_fam
  virt_i = virt[virt$core_domain %in% dom,]
  
  if (any(dom %in% domlist)) {
    # cross tab  
    xtab = table(paste(virt_i$core_domain,virt_i$simarq,sep=" | "), virt_i$tax_factor)
    xtab = xtab[order(rowSums(xtab), decreasing = TRUE),]
    if (is.null(nrow(xtab))) { 
      xtab = t(xtab)
      rownames(xtab) = unique(virt_i$simarq)
    }
    xtab = head(xtab,16)
    barplot(
      t(xtab), 
      col = viridis::plasma(nlevels(virt$tax_factor), end = 0.90, begin = 0.05),
      horiz = TRUE,
      names.arg = paste(str_trunc(rownames(xtab),60), " n=",rowSums(xtab),sep = ""  ),
      las=1, cex.names = 0.6, cex.axis = 0.6, border = NA,ylim = c(0,16), xlim = c(0,160))
    legend(
      "topright", 
      fill = viridis::plasma(nlevels(virt$tax_factor), end = 0.90, begin = 0.05),
      legend = levels(virt$tax_factor), bty="n", cex = 0.2, border = NA)
    title(sprintf("%s\nn=%i",cli, sum(xtab)), cex.main=0.8)
  }
  
}
dev.off()





# closest orthology group per architecture, for euk-like viral homologs
pdf("results_viruses/summary_viral_homologs_tables-arch-ortho.pdf", height = 7, width = 12)
layout(mat=matrix(1:9,nrow=3))
par(mar = c(5.1, 20, 4.1, 2.1))
for (cli in clalist) {
  
  # subset  
  dom = cla[cla$gene_type == cli, ]$gene_fam
  virt_i = virt[virt$core_domain %in% dom,]
  # virt_i = virt_i[virt_i$closest_group == "euk",]
  virt_i$closest_OG_euk_factor = factor(virt_i$closest_OG_euk_filt, levels = unique(virt_i$closest_OG_euk_filt))
  
  if (any(dom %in% domlist)) {
    # cross tab  
    xtab = table(paste(virt_i$core_domain,virt_i$simarq,sep=" | "), virt_i$closest_OG_euk_filt)
    xtab = xtab[order(rowSums(xtab), decreasing = TRUE),levels(virt_i$closest_OG_euk_factor)]
    if (is.null(nrow(xtab))) { 
      xtab = t(xtab)
      rownames(xtab) = unique(virt_i$simarq)
    }
    xtab = head(xtab,16)
    barplot(
      t(xtab), 
      col = rainbow(nlevels(virt_i$closest_OG_euk_factor), end = 0.90),
      horiz = TRUE,
      names.arg = paste(str_trunc(rownames(xtab),60), " n=",rowSums(xtab),sep = ""  ),
      las=1, cex.names = 0.6, cex.axis = 0.6, border = NA,ylim = c(0,16), xlim = c(0,160))
    legend(
      "topright", 
      fill = rainbow(nlevels(virt_i$closest_OG_euk_factor), end = 0.90),
      legend = str_trunc(levels(virt_i$closest_OG_euk_factor),60), bty="n", cex = 0.2, border = NA)
    title(sprintf("%s\nn=%i",cli, sum(xtab)), cex.main=0.8)
  }
  
}
dev.off()

# closest orthology group per architecture, for euk-like viral homologs
pdf("results_viruses/summary_viral_homologs_tables-arch-orthovirfam.pdf", height = 7, width = 12)
layout(mat=matrix(1:9,nrow=3))
par(mar = c(5.1, 20, 4.1, 2.1))
for (cli in clalist) {
  
  # subset  
  dom = cla[cla$gene_type == cli, ]$gene_fam
  virt_i = virt[virt$core_domain %in% dom,]
  virt_i = virt_i[virt_i$closest_group == "euk",]
  virt_i$closest_OG_euk_factor = factor(virt_i$closest_OG_euk_filt, levels = unique(virt_i$closest_OG_euk_filt))
  
  if (any(dom %in% domlist)) {
    # cross tab  
    xtab = table(paste(virt_i$core_domain,virt_i$simarq,sep=" | "), virt_i$tax_factor)
    xtab = xtab[order(rowSums(xtab), decreasing = TRUE),]
    if (is.null(nrow(xtab))) { 
      xtab = t(xtab)
      rownames(xtab) = unique(virt_i$simarq)
    }
    xtab = head(xtab,16)
    barplot(
      t(xtab), 
      col = viridis::plasma(nlevels(virt$tax_factor), end = 0.90, begin = 0.05),
      horiz = TRUE,
      names.arg = paste(str_trunc(rownames(xtab),60), " n=",rowSums(xtab),sep = ""  ),
      las=1, cex.names = 0.6, cex.axis = 0.6, border = NA,ylim = c(0,16), xlim = c(0,160))
    legend(
      "topright", 
      fill = viridis::plasma(nlevels(virt$tax_factor), end = 0.90, begin = 0.05),
      legend = levels(virt$tax_factor), bty="n", cex = 0.2, border = NA)
    title(sprintf("%s\nn=%i",cli, sum(xtab)), cex.main=0.8)
  }
  
}
dev.off()




# diagnostic barplots: how many unique viral genes are classified as 'close to' each other domain, with which confidence?
pdf("results_viruses/summary_viral_homologs_diagnostics.pdf", height = 6, width = 10)
layout(mat=matrix(1:4,nrow=2))
xtab=table(virt$closest_group_top10frac)
barplot(xtab, col="blue", border = "white", ylab = "Freq", 
        main=sprintf("Fraction of hits in closest group, all n=%i", sum(xtab)),
        names.arg = paste(names(xtab), "\nn =",xtab), cex.names = 0.5)

xtab=table(virt[virt$closest_group_factor == "euk",]$closest_group_top10frac)
barplot(xtab, col="blue", border = "white", ylab = "Freq",
        main=sprintf("euk only n=%i", sum(xtab)),
        names.arg = paste(names(xtab), "\nn =",xtab), cex.names = 0.5)

xtab=table(virt[virt$closest_group_factor == "arc",]$closest_group_top10frac)
barplot(xtab, col="blue", border = "white", ylab = "Freq",
        main=sprintf("arc only n=%i", sum(xtab)),
        names.arg = paste(names(xtab), "\nn =",xtab), cex.names = 0.5)

xtab=table(virt[virt$closest_group_factor == "bac",]$closest_group_top10frac)
barplot(xtab, col="blue", border = "white", ylab = "Freq", 
        main=sprintf("bac only n=%i", sum(xtab)),
        names.arg = paste(names(xtab), "\nn =",xtab), cex.names = 0.5)

xtab=table(virt[virt$closest_group_factor == "euk",]$closest_OG_euk_top10frac)
barplot(xtab, col="blue", border = "white", ylab = "Freq", 
        main=sprintf("Fraction of hits to closest orthogroup, euk only n=%i", sum(xtab)),
        names.arg = paste(names(xtab), "\nn =",xtab), cex.names = 0.5)
dev.off()
