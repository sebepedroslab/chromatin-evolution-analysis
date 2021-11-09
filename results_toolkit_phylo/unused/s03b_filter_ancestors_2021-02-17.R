#### Input ####
# load libraries
library(ape)
library(stringr)

pres_fn = "orthogroups_euk.possom.dollo_pres.csv"
summ_fn = "orthogroups_euk.possom.dollo_summary.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
outp_fn = "results_evol"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"

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
mat_summ = merge(mat_summ, cla, by.x = "og_class", by.y = "gene_fam", all.x = T, all.y = F)
rownames(mat_summ) = mat_summ$og_id
mat_summ = mat_summ[mat_summ_original_order,]


# load ancestral presences
mat_pres = read.table(pres_fn, sep="\t", header = T, row.names = 1)
# add num species presence (binarised)
mat_summ$n_presence_bin = apply(mat_pres[,phyl$tip.label], MARGIN = 1, function(x) sum(x>0))
mat_summ$loss_per_presence = mat_summ$n_losses / mat_summ$n_presence_bin


plot(mat_summ$n_presence_bin, mat_summ$n_losses, col="blue", 
     xlab="# species", ylab="# losses")


plot(mat_summ$n_presence_bin, mat_summ$loss_per_presence, col="blue", log="y",
     xlab="# species", ylab="# losses")
abline(h=1, lty=2)

hist(mat_summ$loss_per_presence, breaks = 60, col="blue", border="white")
abline(v=1, lty=2)





stop("ARA")
