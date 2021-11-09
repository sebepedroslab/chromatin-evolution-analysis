#### Input ####
# load libraries
library(ape)
library(stringr)
library(pheatmap)

pres_fn = "orthogroups_euk.ancestral.posteriors_pres.csv"
summ_fn = "orthogroups_euk.ancestral.posteriors_summ.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
outp_fn = "results_evol/"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"

prob_thr = 0.9

graphics.off()

#### Define input ####

# load data
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

gene_counts = table(dat$gene)
genes_repeated = names(gene_counts[gene_counts>1])
ogs_lis_w_dup_genes = dat[dat$gene %in% genes_repeated,"orthogroup"]
dat_dups_full = dat[dat$orthogroup %in% ogs_lis_w_dup_genes,]

# list ogs
ogs_lis = unique(dat_dups_full$orthogroup)
ogs_lis_black = c()

# read species tree
phyl = ape::read.tree(file = phyl_fn)
phyl$edge.length = rep(1, nrow(phyl$edge))
# phyl = phytools::rotateNodes(phyl, nodes = "all")
# dataframe of edges
phyl_edge           = as.data.frame(phyl$edge)
colnames(phyl_edge) = c("edge_start","edge_end")
phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
# dataframe of nodes
phyl_nods = data.frame(taxa = c(phyl$tip.label, phyl$node.label))
phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)

# merge them
phyl_edge = merge(phyl_edge, phyl_nods, all.x = T, all.y = T, by.x = "edge_end", by.y = "edge_end")


# summary of OG data
mat_summ = read.table(summ_fn, sep="\t", header = T)
# remove blacklisted OGs (redundant with others)
mat_summ = mat_summ[!(mat_summ$orthogroup %in% ogs_lis_black),]

# fix OG ids and classes
mat_summ$og_class = stringr::str_split(mat_summ$orthogroup,pattern = "\\.", simplify = T)[,1]
mat_summ$og_id    = stringr::str_split(mat_summ$orthogroup,pattern = ":", simplify = T)[,1]
# fix OG names
mat_summ$og_tmp   = stringr::str_remove(mat_summ$orthogroup, pattern = ":likeclu:.*")
mat_summ$og_tmp   = stringr::str_replace(mat_summ$og_tmp,    pattern = ":like:", replacement = ":like_")
mat_summ$og_name  = stringr::str_split(mat_summ$og_tmp,pattern = ":", simplify = T)[,2]
mat_summ$og_name  = stringr::str_replace(mat_summ$og_name,    pattern = "^like_", replacement = "like:")
mat_summ$og_name_short = stringr::str_trunc(mat_summ$og_name, width = 50)
mat_summ$og_id_short   = paste(mat_summ$og_id, mat_summ$og_name_short, sep=" | ")
rownames(mat_summ) = mat_summ$og_id
mat_summ_original_order = rownames(mat_summ)
# add enzyme classes
# mat_summ = merge(mat_summ, cla, by.x = "og_class", by.y = "gene_fam", all.x = T, all.y = F)
# rownames(mat_summ) = mat_summ$og_id
# mat_summ = mat_summ[mat_summ_original_order,]

# load ancestral presences
mat_pres = read.table(pres_fn, sep="\t", header = T, row.names = 1)
mat_pres_anc = mat_pres [ grep("^SET", rownames(mat_pres)), colnames(mat_pres) %in% phyl$node.label ]

col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

pdf("results_evol/ancestral_presences_table.pdf",width = 20, height = 20)

ogi_alike_heatmap = pheatmap(
  mat_pres_anc,  
  color = col_blue(20), breaks = seq(0,1,length.out = 21)-0.01, 
  cellwidth = 4, cellheight = 4, na_col = "grey", 
  cluster_rows = T, cluster_cols = F, 
  fontsize = 5, 
  legend = F,
  border_color = "white")



dev.off()
