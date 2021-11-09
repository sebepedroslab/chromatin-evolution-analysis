#### Input ####
# load libraries
library(stringr)
library(pheatmap)

dat_fn  = "orthogroups_euk.ancestral.posteriors_pres.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"


#### Define input ####

# load data
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

# read species tree
phyl = ape::read.tree(file = phyl_fn)
phyl$edge.length = rep(1, nrow(phyl$edge))
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



keep_sps1 = tax[tax$Group %in% c("Bilateria"),"Species"]
keep_sps3 = tax[tax$Group %in% c("Fungi"),"Species"]
keep_sps2 = tax[tax$Group %in% c("Streptopoohyta"),"Species"]
keep_anc = c("Metazoa","Fungi","Opisthokonta","Eukaryota")

col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

dat_i = dat[grep("^Chromo", rownames(dat)),]

dat_i = dat_i[,c(keep_sps1,keep_sps2,keep_sps3,keep_anc)]
dat_i = dat_i[rowSums(dat_i[,c(keep_sps1,keep_sps2,keep_sps3)])>0,]
dat_i = dat_i[order(apply(dat_i, 1, which.max), dat_i[ncol(dat_i)],  decreasing = F),]

pdf("unused/matrix_presence.pdf", height = 20, width = 8)
pheatmap(dat_i,
         color = col_blue(10), breaks = seq(0, 1, length.out = 11)-0.01, 
         cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5, gaps_col = c(
             length(keep_sps1),
             length(c(keep_sps1,keep_sps2)),
             length(c(keep_sps1,keep_sps2,keep_sps3))),
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F, number_format = "%i")

dev.off()

