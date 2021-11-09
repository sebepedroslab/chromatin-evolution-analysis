# libraries
library(stringr)
library(pheatmap)
library(ape)
library(igraph)

#### Input ####

# input
clas_fn = "../histone_classification/euk.Histone.to_histdb.Spring-Naive.csv"
alis_fn = "../histone_classification/euk.Histone.to_histdb.csv.gz"
hist_fn = "../../results_toolkit_phylo/gene_counts/euk_genecounts.csv"
phyl_fn = "../../results_toolkit_phylo/data/species_tree.newick"
taxe_fn = "../../data/euk_taxonomy_annotated_2020-08-11.csv"

# load tables
cla = read.table(clas_fn, header = T, sep = "\t", stringsAsFactors = F)
ali = read.table(
    gzfile(alis_fn), sep = "\t", 
    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
    stringsAsFactors = F)
his = read.table(hist_fn, header = F, sep = "\t", stringsAsFactors = F, col.names = c("gene","gene_family"))
his = his[his$gene_family == "Histone",]

# order species (from phylogeny)
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label

# find best alignment
ali = ali[!duplicated(ali$qseqid),]
# ali$histone_class = stringr::str_split(ali$sseqid, pattern = "\\|", simplify = T)[,2]
ali = subset(ali, select = c("qseqid","pident"))

# add classification to histone list
cla$gene = gsub("_\\d+-\\d+$", "", cla$members)
clg = subset(cla, select = c("gene","classification"))
his = merge(his, cla, by.x="gene", by.y= "gene", all.x = T, all.y = F)

# remove histones that can't be classified
his = his[!is.na(his$classification),]
his = his[his$classification != "",]

# add similarity to best hit
his = merge(his, ali, by.x ="members", by.y = "qseqid", all.x = T, all.y = F)


#### Loop: check neighbours ####


# upstream and dowsntream, all together...
nei = data.frame()
list_files_u = list.files(path = "data/", pattern = "*.closest.upstream.bed",full.names = T)
list_files_d = list.files(path = "data/", pattern = "*.closest.downstream.bed",full.names = T)
for (fn in c(list_files_d,list_files_u)) {
    
    # load gene positions
    nei_i = read.table(
        fn, header = F, 
        col.names = c("chr_i","start_i","end_i","transcript_i","strand_i","cluster_i","chr_j","start_j","end_j","transcript_j","strand_j","dist_ij"))
    nei_i = subset(nei_i, select = c("transcript_i","strand_i","transcript_j","strand_j","dist_ij"))
    
    # load old to new gene ID mappings
    dict_getx = read.table(
        sprintf("%s.csv", gsub(pattern = "\\.annot\\.closest\\..*", "",fn)),
        col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
        stringsAsFactors = F)
    dict_getx = dict_getx[!duplicated(dict_getx$sseqid),]
    dict_getx = subset(dict_getx, select = c("qseqid","sseqid"))
    
    # add new ID mappings (histones only!)
    nei_i = merge(nei_i, dict_getx, by.x = "transcript_i", by.y = "sseqid", all.x = T, all.y = F)
    colnames(nei_i)[length(colnames(nei_i))] = "gene_i"
    nei_i = merge(nei_i, dict_getx, by.x = "transcript_j", by.y = "sseqid", all.x = T, all.y = F)
    colnames(nei_i)[length(colnames(nei_i))] = "gene_j"
    
    # concatenate
    nei = rbind(nei, nei_i)
    
}

# discard clusters where there are no neighbour genes
nei = nei[nei$gene_j != ".",]

# which genes are classifiable histones?
nei$is_hist_i = nei$gene_i %in% clg$gene
nei$is_hist_j = nei$gene_j %in% clg$gene

# retain only rows where i and its neighbour are histones
nef = nei[nei$is_hist_i & nei$is_hist_j , ]
nef = subset(nef, select = c("gene_i", "strand_i","gene_j","strand_j","dist_ij"))
nef = merge(nef, clg, by.x = "gene_i", by.y = "gene", all.x = T, all.y = F)
colnames(nef)[length(colnames(nef))] = "class_i"
nef = merge(nef, clg, by.x = "gene_j", by.y = "gene", all.x = T, all.y = F)
colnames(nef)[length(colnames(nef))] = "class_j"

# remove pairs with unclassified members
nef = nef[nef$class_i != "" | nef$class_j != "", ]

# remove pair redundancy
nef$pair = paste(nef$gene_i, nef$gene_j)
nef$pair = unlist(lapply(nef$pair, function(x) paste(sort(unlist(strsplit(x, split = " "))),collapse=",") ) )
nef = nef[!duplicated(nef$pair),]

# pair of histone classes

# first option: relative distance-aware, keep relative order
# nef$pair_hists = NA
# nef[nef$dist_ij >= 0,"pair_hists"] = paste(
#     nef[nef$dist_ij >= 0,]$class_i,
#     nef[nef$dist_ij >= 0,]$class_j,
#     sep = "/")
# nef[nef$dist_ij <  0,"pair_hists"] = paste(
#     nef[nef$dist_ij <  0,]$class_j,
#     nef[nef$dist_ij <  0,]$class_i,
#     sep = "/")
# second option: keep relative order
# nef$pair_hists = paste(nef$class_i,  nef$class_j, sep = "/")
# good option: ignore relative order
nef$pair_hists = paste(nef$class_i, nef$class_j)
nef$pair_hists = unlist(lapply(nef$pair_hists, function(x) paste(sort(unlist(strsplit(x, split = " "))),collapse=",") ) )


# relative orientation of histone pairs
nef$pair_orientation = NA
nef[nef$strand_i == "+" & nef$strand_j == "+" & nef$dist_ij >= 0, "pair_orientation"] = "ht"
nef[nef$strand_i == "+" & nef$strand_j == "-" & nef$dist_ij >= 0, "pair_orientation"] = "tt"
nef[nef$strand_i == "-" & nef$strand_j == "+" & nef$dist_ij >= 0, "pair_orientation"] = "hh"
nef[nef$strand_i == "-" & nef$strand_j == "-" & nef$dist_ij >= 0, "pair_orientation"] = "ht"
nef[nef$strand_i == "+" & nef$strand_j == "+" & nef$dist_ij < 0,  "pair_orientation"] = "ht"
nef[nef$strand_i == "+" & nef$strand_j == "-" & nef$dist_ij < 0,  "pair_orientation"] = "hh"
nef[nef$strand_i == "-" & nef$strand_j == "+" & nef$dist_ij < 0,  "pair_orientation"] = "tt"
nef[nef$strand_i == "-" & nef$strand_j == "-" & nef$dist_ij < 0,  "pair_orientation"] = "ht"

# relative orientation and pair class
nef$pair_hists_ori = paste(nef$pair_hists, nef$pair_orientation)

# nef[nef$dist_ij >= 0,"histone_pair"] = paste(
#     paste(nef[nef$dist_ij >= 0,]$class_i, nef[nef$dist_ij >= 0,]$strand_i, sep = ""), 
#     paste(nef[nef$dist_ij >= 0,]$class_j, nef[nef$dist_ij >= 0,]$strand_j, sep = ""),
#     sep = "/")
# nef[nef$dist_ij <  0,"histone_pair"] = paste(
#     paste(nef[nef$dist_ij <  0,]$class_j, nef[nef$dist_ij <  0,]$strand_j, sep = ""),
#     paste(nef[nef$dist_ij <  0,]$class_i, nef[nef$dist_ij <  0,]$strand_i, sep = ""), 
#     sep = "/")

# add species id
nef$species = stringr::str_split(nef$gene_j , pattern = "_", simplify = T)[,1]
nef$species = factor(nef$species, levels = sps_order)

# remove duplicate entries in each species?
nef_nodups = nef[!duplicated(paste(nef$pair_hists_ori, nef$species)),]


#### Identify clusters ####

# identify connected components using igraph
# (transitive groups from pairs of genes in the table)
nef_s = structure(
    c(nef$gene_i,  nef$gene_j), 
    .Dim = c(nrow(nef), 2L),
    .Dimnames = list(NULL, c("i", "j")))

nef_g = graph.edgelist( as.matrix(nef_s) )
nef_gc = clusters(nef_g)
nef_gcd = data.frame(gene = names(nef_gc$membership), genomic_cluster = nef_gc$membership)

# add genomic cluster info to original dataframe
nef_gcd = merge(nef_gcd, cla, by.x = "gene", by.y = "gene", all.x = T, all.y = F)
nef_gcd$species = stringr::str_split(nef_gcd$gene , pattern = "_", simplify = T)[,1]
nef_gcd$species = factor(nef_gcd$species, levels = sps_order)


# aggregate: list histones in each cluster, and count number of genes
nef_gcd_agg_class = aggregate(nef_gcd$classification, by=list(genomic_cluster = nef_gcd$genomic_cluster,species = nef_gcd$species), FUN = function(x) paste(sort(unique(x)), collapse = ",") )
nef_gcd_agg_count = aggregate(nef_gcd$classification, by=list(genomic_cluster = nef_gcd$genomic_cluster,species = nef_gcd$species), FUN = function(x) length(x) )
nef_gcd_agg_gelis = aggregate(nef_gcd$gene, by=list(genomic_cluster = nef_gcd$genomic_cluster,species = nef_gcd$species), FUN = function(x) paste(x, collapse = ",") )
nef_agg_gen = cbind(nef_gcd_agg_class, nef_gcd_agg_count$x, nef_gcd_agg_gelis$x)
colnames(nef_agg_gen) = c("genomic_cluster","species","members","num_genes","list_genes")
# nef_agg_gen = subset(nef_agg_gen, select = c("species","members","num_genes"))
nef_agg_gen$species = factor(nef_agg_gen$species, levels = sps_order)
# remove singleton clusters
nef_agg_gen = nef_agg_gen[nef_agg_gen$num_genes > 1,]

# aggregate: list types of clusters and count number of genes
nef_agg_clu = aggregate(nef_agg_gen$num_genes, by=list(members = nef_agg_gen$members, species = nef_agg_gen$species), sum)

#### Plots ####

# plot number of histones in each cluster
pdf("cluster_size.pdf", height = 6, width = 5)
barplot(
    table(nef_gc$membership), 
    ylab="num histones in cluster",
    xlab="histone clusters (index)", border = NA)
dev.off()

# plot frequency
pdf("pairs_freq_per_type.pdf", height = 8, width = 5)
par(mar = c(5.1, 16, 4.1, 2.1))
# pair types
tab = sort(table(nef$pair_hists))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (total counts)\nn=%i", sum(tab)))
tab = sort(table(nef_nodups$pair_hists))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (sps-level presence)\nn=%i", sum(tab)))
# orientations
tab = sort(table(nef$pair_orientation))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (total counts)\nn=%i", sum(tab)))
tab = sort(table(nef_nodups$pair_orientation))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (sps-level presence)\nn=%i", sum(tab)))

# types and orientations
tab = sort(table(nef$pair_hists_ori))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (total counts)\nn=%i", sum(tab)), cex.names = 0.6)
tab = sort(table(nef_nodups$pair_hists_ori))
barplot(
    tab, horiz = T, las=1,
    names.arg = paste(names(tab),"n =", tab),
    main = sprintf("histone pairs (sps-level presence)\nn=%i", sum(tab)), cex.names = 0.6)

dev.off()

# plot combinations per species
col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","midnightblue"))
taxe = read.table(taxe_fn, sep = "\t", header = T, stringsAsFactors = T)
taxe = taxe[match(sps_order, taxe$Species),]
sps_gap_ixs = c(1,1+which(diff(as.numeric(taxe$Group))!=0)) - 1

pdf("pairs_freq_per_species.pdf", height = 4, width = 14)
tab = table(nef$pair_hists_ori, nef$species)
tam = matrix(tab, nrow = nrow(tab))
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
tam = tam[order(rowSums(tam), decreasing = T),]
pheatmap(tam,
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         gaps_col = sps_gap_ixs,
         main = sprintf("histone pairs (counts per species)"),
         border_color = "white", display_numbers = T, number_format = "%i")

# without within-species duplicates (only presence/absence)
tab = table(nef_nodups$pair_hists_ori, nef_nodups$species)
tam = matrix(tab, nrow = nrow(tab))
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
tam = tam[order(rowSums(tam), decreasing = T),]
pheatmap(tam,
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         gaps_col = sps_gap_ixs,
         main = sprintf("histone pairs (presence per species)"),
         border_color = "white", display_numbers = T, number_format = "%i")

dev.off()


# plot clusters per species
pdf("cluster_membership_per_species.pdf", height = 4, width = 14)

# num clusters of each type
tab = table(nef_agg_gen$members,nef_agg_gen$species)
tam = matrix(tab, nrow = nrow(tab))
rownames(tam) = rownames(tab)
colnames(tam) = colnames(tab)
order_cluster_type = order(rowSums(tam), decreasing = T)
tam = tam[order_cluster_type,]
pheatmap(tam,
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         gaps_col = sps_gap_ixs,
         main = sprintf("histone cluster membership (num clusters)"),
         border_color = "white", display_numbers = T, number_format = "%i")

# num genes in clusters of each type
tab = aggregate(nef_agg_gen$num_genes, by=list(members = nef_agg_gen$members, species = nef_agg_gen$species), sum)
taw = xtabs(x ~ members + species, data = tab)
tam = matrix(taw, nrow = nrow(taw))
rownames(tam) = rownames(taw)
colnames(tam) = colnames(taw)
tam = tam[order_cluster_type,]
pheatmap(tam,
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         gaps_col = sps_gap_ixs,
         main = sprintf("histone cluster membership (num genes in clusters)"),
         border_color = "white", display_numbers = T, number_format = "%i")


dev.off()



#### Tables ####

# output various summary tables
write.table(nef_agg_gen, file = "summary_clusters.csv", sep="\t", quote = F, row.names = F)
write.table(nef_gcd, file = "summary_clusters_pergene.csv", sep="\t", quote = F, row.names = F)
write.table(nef, file = "summary_pairs.csv", sep="\t", quote = F, row.names = F)
