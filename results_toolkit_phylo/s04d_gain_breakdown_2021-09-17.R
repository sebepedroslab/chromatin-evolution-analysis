#### Input ####
# load libraries
library(ape)
library(stringr)

pres_fn = "orthogroups_euk.ancestral.posteriors_pres.csv"
gain_fn = "orthogroups_euk.ancestral.posteriors_gain.csv"
summ_fn = "orthogroups_euk.ancestral.posteriors_summ.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
outp_fn = "results_evol/"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"

prob_thr = 0.5

graphics.off()

#### Define input ####

# load data
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE)
cla = read.table(clas_fn, header = FALSE, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_counts = table(dat$gene)
genes_repeated = names(gene_counts[gene_counts > 1])
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
phyl_edge = merge(phyl_edge, phyl_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")


# summary of OG data
mat_summ = read.table(summ_fn, sep="\t", header = TRUE)
# remove blacklisted OGs (redundant with others)
mat_summ = mat_summ[!(mat_summ$orthogroup %in% ogs_lis_black),]

# fix OG ids and classes
mat_summ$og_class = stringr::str_split(mat_summ$orthogroup,pattern = "\\.", simplify = TRUE)[,1]
mat_summ$og_id    = stringr::str_split(mat_summ$orthogroup,pattern = ":", simplify = TRUE)[,1]
# fix OG names
mat_summ$og_tmp   = stringr::str_remove(mat_summ$orthogroup, pattern = ":likeclu:.*")
mat_summ$og_tmp   = stringr::str_replace(mat_summ$og_tmp,    pattern = ":like:", replacement = ":like_")
mat_summ$og_name  = stringr::str_split(mat_summ$og_tmp,pattern = ":", simplify = TRUE)[,2]
mat_summ$og_name  = stringr::str_replace(mat_summ$og_name,    pattern = "^like_", replacement = "like:")
mat_summ$og_name_short = stringr::str_trunc(mat_summ$og_name, width = 50)
mat_summ$og_id_short   = paste(mat_summ$og_id, mat_summ$og_name_short, sep=" | ")
rownames(mat_summ) = mat_summ$og_id
mat_summ_original_order = rownames(mat_summ)

# load ancestral presences
mat_pres = read.table(pres_fn, sep="\t", header = TRUE, row.names = 1)
mat_gain = read.table(gain_fn, sep="\t", header = TRUE, row.names = 1)

# add num species presence (binarised)
mat_summ$n_presence_bin = apply(mat_pres[,phyl$tip.label], MARGIN = 1, function(x) sum(x > 0))

# simplify OG names
rownames(mat_pres) = mat_summ$og_id
rownames(mat_gain) = mat_summ$og_id

# remove blacklisted OGs
print(sprintf("Remove %i OGs", length(ogs_lis_black)))
mat_summ = mat_summ[!(rownames(mat_summ) %in% ogs_lis_black),]
mat_pres = mat_pres[!(rownames(mat_pres) %in% ogs_lis_black),]
mat_gain = mat_gain[!(rownames(mat_gain) %in% ogs_lis_black),]

# order taxa according to phylogeny?
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]

# create factors for taxa and families
# mat_summ$gain = factor(mat_summ$gain, levels = taxa_lis)
mat_summ$og_class = factor(mat_summ$og_class, levels = unique(mat_summ$og_class))

# add classes to per-gene table
dat = merge(dat, mat_summ, by.x="orthogroup", by.y="orthogroup", all.x = TRUE, all.y = FALSE)
# add species 
dat$taxa = stringr::str_split(dat$gene,pattern = "_", simplify = TRUE)[,1]
dat$taxa = factor(dat$taxa, levels = taxa_lis)
# remove blacklisted ogs
dat = dat[!(dat$og_class %in% ogs_lis_black), ]

# list of gene types
gene_types = unique(as.character(mat_summ$gene_type))
gene_classes = unique(as.character(mat_summ$gene_class))

# list of nodes of interest
list_nodes = c("Metazoa","Holozoa","Bilateria","BilCniTri","Opisthokonta","Amorphea","Diaphoratickes","Fungi")
list_nodes = phyl$node.label

# na to 0
mat_gain [ is.na(mat_gain) ] = 0



#### Node-specific gains ####

# loop through gene types
for (type in gene_types) {
  
  message(sprintf("gains per node | %s", type))
  
  pdf(sprintf("results_evol/breakdown_gains_per_node.type.%s.pdf", type),width = 6.5, height = 6.5)
  layout(matrix(1:16, nrow = 4))
  
  for (node in list_nodes) {
    
    # restrict analysis to OGs present in LMCA
    mat_summ_i = mat_summ[mat_gain[,node] > prob_thr ,]
    mat_summ_i$og_class_factor = factor(mat_summ_i$og_class, levels = as.character(cla$gene_fam))
    
    # table    
    support_summary = table(mat_summ_i[mat_summ_i$gene_type == type,]$og_class)
    support_labels = sprintf("%s n=%i", names(support_summary),support_summary)
    support_labels [ grepl("n=0$", support_labels)] = ""
    if (sum(support_summary) > 0) {
      pie(support_summary, las=1, border = "white", lwd = 0.6,
          col = rainbow(nlevels(mat_summ_i$og_class_factor)),
          cex=0.6, main=sprintf("%s gain in %s",type, node), cex.main=0.8,cex.sub=0.8,
          labels = support_labels,
          sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
    }
    
  }
  # legend
  pie(1, col = NA, border = NA, labels=NA)
  legend("topright", legend=levels(mat_summ_i$og_class), fill=rainbow(nlevels(mat_summ_i$og_class_factor)), cex=0.3, bty = "n")
  dev.off()
  
}


# loop through gene classes
for (type in gene_classes) {
  
  message(sprintf("gains per node | %s", type))
  
  pdf(sprintf("results_evol/breakdown_gains_per_node.class.%s.pdf", type),width = 6.5, height = 6.5)
  layout(matrix(1:16, nrow = 4))
  
  for (node in list_nodes) {
    
    # restrict analysis to OGs present in LMCA
    mat_summ_i = mat_summ[mat_gain[,node] > prob_thr ,]
    mat_summ_i$og_class_factor = factor(mat_summ_i$og_class, levels = as.character(cla$gene_fam))
    
    # table    
    support_summary = table(mat_summ_i[mat_summ_i$gene_class == type,]$og_class)
    support_labels = sprintf("%s n=%i", names(support_summary),support_summary)
    support_labels [ grepl("n=0$", support_labels)] = ""
    if (sum(support_summary) > 0) {
      pie(support_summary, las=1, border = "white", lwd = 0.6,
          col = rainbow(nlevels(mat_summ_i$og_class_factor)),
          cex=0.6, main=sprintf("%s gain in %s",type, node), cex.main=0.8,cex.sub=0.8,
          labels = support_labels,
          sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
    }
    
  }
  # legend
  pie(1, col = NA, border = NA, labels=NA)
  legend("topright", legend=levels(mat_summ_i$og_class), fill=rainbow(nlevels(mat_summ_i$og_class_factor)), cex=0.3, bty = "n")
  dev.off()
  
}




#### Lineage-specific gains (node+descendants) ####

# loop through gene types
for (type in gene_types) {
  
  message(sprintf("gains per lineage | %s", type))
  
  pdf(sprintf("results_evol/breakdown_gains_per_lineage.type.%s.pdf", type),width = 6.5, height = 6.5)
  layout(matrix(1:16, nrow = 4))
  
  for (node in list_nodes) {
    
    node_descs = c(phyl$tip.label, phyl$node.label) [ getDescendants(phyl, length(phyl$tip.label) + which(phyl$node.label == node)) ]
    node_descs = c(node, node_descs)
    # restrict analysis to OGs present in LMCA
    ogs_bool = apply(mat_gain[,node_descs] > prob_thr, 1, any) # this will take any OG which has >50% gain probability in at least one of the nodes of interest
    # ogs_bool = rowSums(mat_gain[,node_descs]) > 0.9            # this will take any OG which has >90% gain probability aggregated over all nodes of interest
    mat_summ_i = mat_summ[ ogs_bool ,]
    mat_summ_i$og_class_factor = factor(mat_summ_i$og_class, levels = as.character(cla$gene_fam))
    
    # table    
    support_summary = table(mat_summ_i[mat_summ_i$gene_type == type,]$og_class)
    support_labels = sprintf("%s n=%i", names(support_summary),support_summary)
    support_labels [ grepl("n=0$", support_labels)] = ""
    if (sum(support_summary) > 0) {
      pie(support_summary, las=1, border = "white", lwd = 0.6,
          col = rainbow(nlevels(mat_summ_i$og_class_factor)),
          cex=0.6, main=sprintf("%s gain from %s on",type, node), cex.main=0.8,cex.sub=0.8,
          labels = support_labels,
          sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
    }
    
  }
  # legend
  pie(1, col = NA, border = NA, labels=NA)
  legend("topright", legend=levels(mat_summ_i$og_class), fill=rainbow(nlevels(mat_summ_i$og_class_factor)), cex=0.3, bty = "n")
  dev.off()
  
}
