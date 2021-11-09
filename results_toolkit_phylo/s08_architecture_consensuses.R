# load libraries
library(scales)
library(stringr)
library(pheatmap)
library(ape)
library(tidyr)

graphics.off()

#### Define input ####

clas_fn = "../data/gene_families_hmm.csv"
arqs_fo = "results_domains/"
phyl_fn = "data/species_tree.newick"
tedo_fn = "../data/transposable_element_domains.csv"
tree_fo = "gene_trees//"
garq_fn = "architectures.csv"
gcla_fn = "gene_counts/euk_genecounts.csv"
taxe_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"   # taxonomic info for eukaryotes
dat_fn = "orthogroups_euk.csv"

col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

# load tables
cla = read.table(clas_fn, header = F, sep = "\t", 
                 col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
                 stringsAsFactors = F)

# load orthogroups
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
dat$og_gene_fam = stringr::str_split(dat$orthogroup, "\\.", simplify = T)[,1]

# species tree for species ordering
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label

# load taxonomy and find out where to add spacing in heatmaps
taxe = read.table(taxe_fn, sep = "\t", header = T, stringsAsFactors = T)
taxe = taxe[match(sps_order, taxe$Species),]
sps_gap_ixs = c(1,1+which(diff(as.numeric(taxe$Group))!=0)) - 1

# load architectures
garq = read.table(garq_fn, header = F, sep="\t", col.names = c("gene","architecture"))
gcla = read.table(gcla_fn, header = F, sep = "\t", col.names = c("gene","gene_fam"))
garq_cla = merge(garq, gcla, by.x = "gene", by.y="gene", all.x = F, all.y = T)
garq_cla = merge(garq_cla, cla, by.x = "gene_fam", by.y = "gene_fam", all.x = T, all.y = F)

# simplify architectures (collapse **consecutive** repeats)
# simple_arqs_list = stringr::str_split(garq_cla$architecture, pattern = " ")
# simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
# garq_cla$simarq = simple_arqs_vect

# super-simplify architectures (collaps consecutive repeats, ignore "_num" suffixes)
simple_arqs_list = stringr::str_split(garq_cla$architecture, pattern = " ")
simple_arqs_list = sapply(simple_arqs_list, function(x) gsub(pattern = "_\\d+$", "", x))
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
garq_cla$simarq = simple_arqs_vect


# add species
garq_cla$species = stringr::str_split(garq_cla$gene, pattern = "_", simplify = T)[,1]

list_fams = unique(cla$gene_type)
list_fams = c("Acetylase","Deacetylase","Methylase","Demethylase","Readers","Remodeller","Chaperones")
# list_fams = c("Methylase")

for (fam in list_fams) {
  
  # subset architectures to gene type
  garq_i = garq_cla[garq_cla$gene_type == fam,]
  fams_in_type = unique(cla[cla$gene_type == fam,"gene_fam"])
  
  # add orthogroup
  garq_i = merge(garq_i, dat, by.x = "gene", by.y = "gene", all.x = T, all.y = F)
  
  # exclude orthogroups from genefams not in this gene type
  garq_i = garq_i [ garq_i$og_gene_fam %in% fams_in_type, ]
  
  # subset OGs of interest
  garq_i = garq_i [ !is.na(garq_i$orthogroup) , ]
  garq_i = garq_i [ !grepl(":like:",garq_i$orthogroup) , ]
  garq_i = garq_i [ !grepl(":NA$",garq_i$orthogroup) , ]
  
  # list of OGs
  coun_ogs = sort(table(garq_i$orthogroup), decreasing = T)
  list_ogs = names(coun_ogs)
  
  #### Common architectures ####
  # identify which architectures are common (at least 10%) within an orthogroup
  # and plot them uglily, for reference
  yi=0
  arq_counts_ids_list = c()
  pdf(file = sprintf("%s/consenus_architectures_%s.dots.pdf", arqs_fo, fam),width = 16, height = 10)
  plot(NA,NA,xlim=c(0,6), ylim=c(0,40),main=fam, axes=F, xlab=NA, ylab=NA)
  for(ogi in list_ogs) {
    
    yi=yi+1
    
    garq_o = garq_i[garq_i$orthogroup == ogi,]
    garq_u = garq_o[!duplicated(paste(garq_o$species, garq_o$simarq)),]
    num_species = length(unique(garq_o$species))
    num_genes = nrow(garq_o)
    arq_counts = table(garq_o$simarq)
    arq_counts_top = arq_counts[ arq_counts / num_genes > 0.1 ]
    arq_counts_ids = names(arq_counts_top)
    arq_counts_sps = unlist(lapply(arq_counts_ids, function(x) nrow(garq_u[garq_u$simarq == x,] )))
    names(arq_counts_sps) = arq_counts_ids
    arq_counts_spsfrac = arq_counts_sps / num_species
    arq_counts_spsfrac = sort(arq_counts_spsfrac, decreasing = T)
    
    # plot as dotplot
    points(x=1:length(arq_counts_spsfrac), y=rep(yi,length(arq_counts_spsfrac)), cex=sqrt(arq_counts_spsfrac*4), col=alpha("springgreen4",0.5), pch=19)
    text(x=1:length(arq_counts_spsfrac), y=rep(yi,length(arq_counts_spsfrac)), labels = names(arq_counts_spsfrac), adj = 0, col=alpha("darkgreen",0.8),cex=0.6)
    
    arq_counts_ids_list = unique(c(arq_counts_ids_list, arq_counts_ids))
    
  }
  # og names in dotplot
  text(x=0.8, y=1:length(list_ogs), labels = list_ogs, adj=1,cex=0.6)
  dev.off()
  
  # list of common architectures within a given gene family
  arq_counts_ids_list = arq_counts_ids_list[arq_counts_ids_list != ""]
  
  
  # HEATMAP
  # crosstabulate OG and simarq
  xtab = xtabs(formula = ~ orthogroup + simarq, data=garq_i)
  total_genes_per_og = rowSums(xtab)
  xtab = xtab[,arq_counts_ids_list]
  
  # ordering of architectures (following presence in ogs)
  ord_col = order(apply(xtab, 2, which.max), decreasing = T)
  xtab = xtab[,ord_col]
  
  # turn into fraction of archs within OG (total number of genes, not just
  # genes with common architectures)
  xtab = xtab / total_genes_per_og
  
  # standardise rownames and colnames
  xtab_rows = stringr::str_pad(stringr::str_trunc(rownames(xtab), 40, "right"), width = 45, "right")
  xtab_cols = stringr::str_pad(colnames(xtab), 45, "right")
  pdf(file = sprintf("%s/consenus_architectures_%s.table.pdf", arqs_fo, fam),width = 8, height = 7)
  pheatmap(xtab,
           color = col_blue(10), breaks = seq(0,1,length.out = 11),
           cellwidth = 3, cellheight = 3, na_col = "grey",
           # cluster_rows = T, cluster_cols = T,
           # clustering_distance_cols="correlation",
           # clustering_distance_rows="correlation",
           cluster_rows = F, cluster_cols = F,
           fontsize = 4,
           labels_row = xtab_rows,
           labels_col = xtab_cols,
           main = sprintf("%s", fam),
           border_color = "white", display_numbers = F, number_format = "%.2f")
  dev.off()
  
  
  
  #### Common single domains ####
  # identify which domains are common (at least 10%) within an orthogroup
  arq_counts_ids_list = c()
  for(ogi in list_ogs) {
    
    garq_o = garq_i[garq_i$orthogroup == ogi,]
    garq_u = garq_o[!duplicated(paste(garq_o$species, garq_o$simarq)),]
    num_genes = nrow(garq_o)
    domain_vector = stringr::str_split(garq_o$simarq, pattern = ",")
    domain_counts = table(unlist(domain_vector))
    arq_counts_top = domain_counts[ domain_counts / num_genes > 0.1 ]
    arq_counts_ids = names(arq_counts_top)
    arq_counts_ids_list = unique(c(arq_counts_ids_list, arq_counts_ids))
    
  }
  
  # list of common architectures within a given gene family
  arq_counts_ids_list = arq_counts_ids_list[arq_counts_ids_list != ""]
  
  
  # HEATMAP
  # get extended table with one domain per row
  dom_per_gene = stringr::str_split(garq_i$simarq, pattern = ",")
  names(dom_per_gene) = garq_i$gene
  dom_per_gene_table = data.frame(domain = unlist(dom_per_gene))
  dom_per_gene_table$gene = 
    unlist(lapply(
      names(dom_per_gene), 
      function(x) 
        rep(x, length(dom_per_gene[[x]]))
    ))
  dom_per_gene_table = merge(dom_per_gene_table, dat, by.x = "gene", by.y = "gene", all.x = T, all.y = F)
  # exclude orthogroups from genefams not in this gene type
  dom_per_gene_table = dom_per_gene_table [ dom_per_gene_table$og_gene_fam %in% fams_in_type, ]
  # subset OGs of interest
  dom_per_gene_table = dom_per_gene_table [ !is.na(dom_per_gene_table$orthogroup) , ]
  dom_per_gene_table = dom_per_gene_table [ !grepl(":like:",dom_per_gene_table$orthogroup) , ]
  dom_per_gene_table = dom_per_gene_table [ !grepl(":NA$",dom_per_gene_table$orthogroup) , ]
  
  
  # crosstabulate OG and domains
  xtab = xtabs(formula = ~ orthogroup + domain, data=dom_per_gene_table)
  xtab = xtab[,arq_counts_ids_list]
  
  # ordering of architectures (following presence in ogs)
  ord_col = order(apply(xtab, 2, which.max),colSums(xtab),  decreasing = T)
  xtab = xtab[,ord_col]
  
  # turn into fraction of archs within OG (total number of genes, not just
  # genes with common architectures)
  xtab = xtab / total_genes_per_og
  
  # standardise rownames and colnames
  xtab_rows = stringr::str_pad(stringr::str_trunc(rownames(xtab), 40, "right"), width = 45, "right")
  xtab_cols = stringr::str_pad(colnames(xtab), 45, "right")
  gaps_row = c(1,1+which(diff( as.numeric(as.factor(stringr::str_split(xtab_rows, "\\.", simplify = T)[,1])) )!=0)) - 1
  pdf(file = sprintf("%s/consenus_domains_%s.table.pdf", arqs_fo, fam),width = 8, height = 7)
  pheatmap(xtab,
           color = col_blue(10), breaks = seq(0,1,length.out = 11),
           cellwidth = 3, cellheight = 3, na_col = "grey",
           cluster_rows = F, cluster_cols = F,
           # clustering_distance_cols="correlation",
           # clustering_distance_rows="correlation",
           # cluster_rows = F, cluster_cols = F,
           fontsize = 4,
           labels_row = xtab_rows,
           labels_col = xtab_cols,
           # gaps_row = gaps_row,
           main = sprintf("%s", fam),
           border_color = "white", display_numbers = F, number_format = "%.2f")
  
  
  # now add barcodes with species presence
  dom_per_gene_table$species = stringr::str_split(dom_per_gene_table$gene, "_", simplify = T)[,1]
  dom_per_gene_table$species = factor(dom_per_gene_table$species, levels = taxe$Species)
  xtax = xtabs(formula = ~ orthogroup + species, data=dom_per_gene_table)
  pheatmap(xtax,
           color = col_blue(10), breaks = seq(0,1,length.out = 11),
           cellwidth = 0.25, cellheight = 3, na_col = "grey",
           cluster_rows = F, cluster_cols = F,
           # clustering_distance_cols="correlation",
           # clustering_distance_rows="correlation",
           # cluster_rows = F, cluster_cols = F,
           fontsize = 4,
           labels_row = xtab_rows,
           labels_col = NA,
           # gaps_row = gaps_row,
           main = sprintf("%s", fam),
           border_color = NA, display_numbers = F, number_format = "%.2f")
  
  # now same but at the macro level
  dom_per_gene_table_tax = merge(dom_per_gene_table, taxe, by.x = "species", by.y = "Species", all.x = T, all.y = F, sort = F)
  dom_per_gene_table_tax$Group = factor(dom_per_gene_table_tax$Group, levels = unique(as.character(taxe$Group)))
  xtax = xtabs(formula = ~ orthogroup + Group, data=dom_per_gene_table_tax)
  pheatmap(xtax,
           color = col_blue(10), breaks = seq(0,1,length.out = 11),
           cellwidth = 0.75, cellheight = 3, na_col = "grey",
           cluster_rows = F, cluster_cols = F,
           # clustering_distance_cols="correlation",
           # clustering_distance_rows="correlation",
           # cluster_rows = F, cluster_cols = F,
           fontsize = 4,
           labels_row = xtab_rows,
           # gaps_row = gaps_row,
           main = sprintf("%s", fam),
           border_color = NA, display_numbers = F, number_format = "%.2f")
  
  
  dev.off()
  
}





#### Dedicated readers analysis ####

source("../scripts/venn_diagrams.R")

### Get readers
fam="Readers"
# subset architectures to gene type
garq_i = garq_cla[garq_cla$gene_type == fam,]

### venn diagrams of dual classification genes
pdf(file = sprintf("%s/venn_readers_v_catalytic.pdf", arqs_fo),width = 2.5, height = 2.5)

# first, catalytic
nonfam = c("Acetylase","Deacetylase","Methylase","Demethylase","Remodeller")
garq_j = garq_cla[garq_cla$gene_type %in% nonfam,]
ovs = venn.two(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  catname1 = "Readers", catname2 = "Catalytic", 
  col1 = "purple", col2 = "chartreuse3",
  main="Catalytic")

# second, non-readers (joint)
nonfam = unique(cla$gene_type)
nonfam = nonfam [ nonfam != "Readers" ]
garq_j = garq_cla[garq_cla$gene_type %in% nonfam,]
ovs = venn.two(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  catname1 = "Readers", catname2 = "All", 
  col1 = "purple", col2 = "chartreuse3",
  main="All")


# third, one by one
for (faj in nonfam) {
  
  garq_j = garq_cla[garq_cla$gene_type == faj,]
  
  ovs = venn.two(
    unique(garq_i$gene), 
    unique(garq_j$gene), 
    catname1 = "Readers", catname2 = faj, 
    col1 = "purple", col2 = "chartreuse3",
    main=faj)
  
  
}
dev.off()


### opposite: catalytic v readers
pdf(file = sprintf("%s/venn_catalytic_v_readers.pdf", arqs_fo),width = 2.5, height = 2.5)

nonfam = c("Acetylase","Deacetylase","Methylase","Demethylase","Remodeller")
garq_i = garq_cla[garq_cla$gene_type  %in% nonfam, ]
for (faj in unique(garq_cla[garq_cla$gene_type == "Readers",]$gene_fam)) {
  
  garq_j = garq_cla[garq_cla$gene_fam == faj,]
  
  ovs = venn.two(
    unique(garq_i$gene), 
    unique(garq_j$gene), 
    catname1 = "Catalytic", catname2 = faj, 
    col1 = "purple", col2 = "chartreuse3",
    main=faj)
  
  
}
dev.off()



# third, one by one
pdf(file = sprintf("%s/venn_readers_internal.pdf", arqs_fo),width = 2.5, height = 2.5)


# main readers
garq_i = garq_cla[garq_cla$gene_fam == "Chromo",]
garq_j = garq_cla[garq_cla$gene_fam == "PHD",]
garq_k = garq_cla[garq_cla$gene_fam == "Bromodomain",]
ovs = venn.three(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  unique(garq_k$gene), 
  catname1 = "Chromo",
  catname2 = "PHD", 
  catname3 = "Bromo",
  col1 = "purple", col2 = "chartreuse3", col3 = "orange",
  eulerbool = TRUE,
  main = "main")

# chromo v others
garq_i = garq_cla[garq_cla$gene_fam == "Chromo",]
garq_j = garq_cla[!garq_cla$gene_fam %in% c("Chromo","Bromodomain","PHD") & garq_cla$gene_class == "Readers",]
ovs = venn.two(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  catname1 = "Chromo",
  catname2 = "Other", 
  col1 = "purple", col2 = "chartreuse3",
  eulerbool = TRUE,
  main = "main")

# chromo v others
garq_i = garq_cla[garq_cla$gene_fam == "Bromodomain",]
garq_j = garq_cla[!garq_cla$gene_fam %in% c("Chromo","Bromodomain","PHD") & garq_cla$gene_class == "Readers",]
ovs = venn.two(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  catname1 = "Bromodomain",
  catname2 = "Other", 
  col1 = "purple", col2 = "chartreuse3",
  eulerbool = TRUE,
  main = "main")

# chromo v others
garq_i = garq_cla[garq_cla$gene_fam == "PHD",]
garq_j = garq_cla[!garq_cla$gene_fam %in% c("Chromo","Bromodomain","PHD") & garq_cla$gene_class == "Readers",]
ovs = venn.two(
  unique(garq_i$gene), 
  unique(garq_j$gene), 
  catname1 = "PHD",
  catname2 = "Other", 
  col1 = "purple", col2 = "chartreuse3",
  eulerbool = TRUE,
  main = "main")

dev.off()


# hh = hist(table(garq_cla$gene), breaks = 0:5)
