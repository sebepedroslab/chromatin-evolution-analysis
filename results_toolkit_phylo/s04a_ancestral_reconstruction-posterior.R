# load libs
library(stringr)
library(ape)
library(ggtree)
library(scales)
library(pheatmap)
library(phytools)

#### Define input ####


# input
prob_thr = 0.5
dat_fn = "orthogroups_euk.csv"
ort_fn = "orthogroups_euk.ancestral.posteriors.csv"
phy_fn = "data/species_tree.newick"
cla_fn = "../data/gene_families_hmm.csv"
out_fo = "results_evol"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
arq_fn = "architectures.csv"


# read gene classification, taxonomy, etc
cla = read.table(cla_fn, header = FALSE, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE)
arq = read.table(arq_fn, sep = "\t", stringsAsFactors = FALSE, header = FALSE, col.names = c("gene","architecture"))


# read species tree
phyl = ape::read.tree(phy_fn)
sps_list = phyl$tip.label
phyl$edge.length = rep(1, nrow(phyl$edge))
# prepare species dataframe:
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

# heatmap col
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

#### Read Posteriors ####

# read orthogroup data
ort = read.table(ort_fn, header = TRUE)
# drop ASBENT entry (useless here)
ort = ort[ort$Family != "ABSENT", ]
ort_taxa = str_remove(colnames(ort)[grepl("\\.1$", colnames(ort))], "\\.1")
rownames(ort) = ort$Family

# get presence, gains and loss matrices
mat_pre1 = ort[,grepl("\\.1$", colnames(ort))]
mat_prem = ort[,grepl("\\.m$", colnames(ort))]
mat_gain = ort[,grepl("\\.gain$", colnames(ort))]
mat_loss = ort[,grepl("\\.loss$", colnames(ort))]
# NaN to zero
# mat_pre1[is.na(mat_pre1)] = 0
# mat_prem[is.na(mat_prem)] = 0
# mat_loss[is.na(mat_loss)] = 0
# mat_gain[is.na(mat_gain)] = 0
# sum one/more presence probs
mat_pres = mat_pre1 + mat_prem
# add empty ANCESTRAL node to gains and losses
mat_gain = cbind(mat_gain, NA)
mat_loss = cbind(mat_loss, NA)
# colnames
colnames(mat_pres) = ort_taxa
colnames(mat_gain) = ort_taxa
colnames(mat_loss) = ort_taxa
# reorder to match species tree
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]
mat_loss = mat_loss[,taxa_lis]

# get summary matrix
mat_summ = ort[,c("Family","Pattern")]
colnames(mat_summ) = c("orthogroup","Pattern")

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
# add enzyme classes
mat_summ = merge(mat_summ, cla, by.x = "og_class", by.y = "gene_fam", all.x = TRUE, all.y = FALSE)
rownames(mat_summ) = mat_summ$og_id
mat_summ = mat_summ[mat_summ_original_order,]

# drop unwanted (those not in the master gene classification table `cla`)
ixs_keep = which(!is.na(mat_summ$gene_type))
mat_summ = mat_summ [ ixs_keep, ]
mat_pres = mat_pres [ ixs_keep, ]
mat_loss = mat_loss [ ixs_keep, ]
mat_gain = mat_gain [ ixs_keep, ]

# simplify OG names
rownames(mat_gain) = mat_summ$og_id
rownames(mat_loss) = mat_summ$og_id
rownames(mat_pres) = mat_summ$og_id

# add summary info
mat_summ$sum_gain = apply(mat_gain, 1, sum)
mat_summ$sum_loss = apply(mat_loss, 1, sum)
mat_summ$sum_pres_total = apply(mat_pres, 1, sum)
mat_summ$sum_pres_extant = apply(mat_pres[,colnames(mat_pres) %in% sps_list], 1, sum)

# create factors for families
mat_summ$og_class = factor(mat_summ$og_class, levels = unique(mat_summ$og_class))

# add classes to per-gene table
dat = merge(dat, mat_summ, by.x="orthogroup", by.y="orthogroup", all.x = TRUE, all.y = FALSE)
# add species 
dat$taxa = stringr::str_split(dat$gene,pattern = "_", simplify = TRUE)[,1]
dat$taxa = factor(dat$taxa, levels = taxa_lis)

# list of gene types
gene_types = unique(as.character(cla$gene_type))
gene_types = gene_types[gene_types %in% as.character(mat_summ$gene_type)]
gene_classes = unique(as.character(cla$gene_class))
gene_classes = gene_classes[gene_classes %in% as.character(mat_summ$gene_class)]

# fix architecture names
simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
arq$simarq = simple_arqs_vect

##### Output matrices #####

# create presence, gain and loss matrices per node
write.table(mat_pres, "orthogroups_euk.ancestral.posteriors_pres.csv", sep = "\t", quote = FALSE)
write.table(mat_loss, "orthogroups_euk.ancestral.posteriors_loss.csv", sep = "\t", quote = FALSE)
write.table(mat_gain, "orthogroups_euk.ancestral.posteriors_gain.csv", sep = "\t", quote = FALSE)
write.table(mat_summ, "orthogroups_euk.ancestral.posteriors_summ.csv", sep = "\t", quote = FALSE)



##### plot global cladograms: all species #####

# first, plot model phylogenies
pdf(file = sprintf("%s/phylo_global_sps.00.pdf", out_fo),width = 8, height = 8)
ggphy = ggtree(phyl, ladderize = FALSE, layout = "rectangular",branch.length = "none", lwd=0.5, color="gray") +
  geom_nodelab(color="red", size=1.5) +
  ggplot2::scale_y_reverse()
print(ggphy + geom_tiplab(offset=0, size=1.5))
open_tree(ggphy, 5) + geom_tiplab(offset=0.5, size=1.5)
open_tree(ggphy, 180) + geom_tiplab(offset=0.5, size=1.5)
dev.off()



# now aggregated cladograms
for (type in gene_types) {
  
  pdf(file = sprintf("%s/phylo_global-type_sps.%s.pdf", out_fo, type),width = 8, height = 20)
  print(sprintf("Plot aggregated counts %s | species", type))
  
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  #   mat_gpl_sum$gain = colSums(mat_gain[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$loss = colSums(mat_loss[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$pres = colSums(mat_pres[cla_boo,], na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyl_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%#i")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%#i")," / -",sprintf(phyl_data$loss, fmt = "%#i"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then gains
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  edgelabels(text=sprintf("+%i",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-20,20))
  title(sprintf("Gains and losses in ancestral nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot gain and losses barplots: extant  
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[phyl_data$is_tip,]$gain,
      loss = phyl_data[phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
}





##### plot global cladograms: euk taxonomic groups, per gene type #####


# prune species tree until groups of interest are tips:
tax_groups = unique(tax$Group)
phyg = phyl
for (taxi in tax_groups) {
  sps_in_tax_group = tax[tax$Group == taxi, "Species"]
  # drop tips one by one, except for first tip
  tip_to_keep = sps_in_tax_group[1]
  tip_to_drop = sps_in_tax_group[-1]
  for (spi in tip_to_drop) {
    phyg = ape::drop.tip(phy = phyg, tip=spi, trim.internal = TRUE)
  }
  # if group contains more than one species, reassign tip name to 
  # its LCA name. Otherwise, keep tip name (e.g. Ttra = Apusozoa)
  if (length(sps_in_tax_group) > 1) {
    phyg$tip.label[phyg$tip.label == tip_to_keep] = taxi
  }
}
# reassing edge lengths (all 1)
phyg$edge.length = rep(1, length(phyg$edge.length))
ape::write.tree(phyg, file = "data/species_tree_pruned.newick")

# dataframe of edges
phyg_edge           = as.data.frame(phyg$edge)
colnames(phyg_edge) = c("edge_start","edge_end")
phyg_edge$ix_edges = as.numeric(rownames(phyg_edge))
phyg_edge$ends_in_tip = phyg_edge$edge_end <= length(phyg$tip.label)
# dataframe of nodes
phyg_nods = data.frame(taxa = c(phyg$tip.label, phyg$node.label))
phyg_nods$edge_end = as.numeric(rownames(phyg_nods))
phyg_nods$is_tip   = phyg_nods$edge_end <= length(phyg$tip.label)
# merge them
phyg_edge = merge(phyg_edge, phyg_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")


# first, plot model phylogenies
pdf(file = sprintf("%s/phylo_global_groups.00.pdf", out_fo),width = 3, height = 4)
ggphy = ggtree(phyg, ladderize = FALSE, layout = "rectangular",branch.length = "none", lwd=0.5, color="gray") +
  geom_nodelab(color="red", size=1.5) +
  ggplot2::scale_y_reverse()
print(ggphy + geom_tiplab(offset=0, size=1.5))
open_tree(ggphy, 5) + geom_tiplab(offset=0.5, size=1.5)
open_tree(ggphy, 180) + geom_tiplab(offset=0.5, size=1.5)
write.tree(phyg, sprintf("%s/phylo_global_groups.00.newick", out_fo))
dev.off()

# now aggregated cladograms
for (type in gene_types) {
  
  pdf(file = sprintf("%s/phylo_global-type_groups.%s.pdf", out_fo, type),width = 2, height = 3)
  par(mar = c(2, 2, 4.1, 2))
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  #   mat_gpl_sum$gain = colSums(mat_gain[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$loss = colSums(mat_loss[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$pres = colSums(mat_pres[cla_boo,], na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE, all.y=FALSE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%#i")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%#i")," / -",sprintf(phyl_data$loss, fmt = "%#i"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  edgelabels(text=sprintf("+%i",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  dev.off()
  
  
  
  
  # plot distribution of OG presence counts in the descendant tips from each "tip" taxonomic group
  tip_distributions = list()
  tip_counts = list()
  for (tip in phyg$tip.label) {
    
    # if display tip is an extant species (will happen for single-species taxon groups):
    if (tip %in% phyl$tip.label) {
      tip_descendants = tip
    } else {
      tip_descendants = tax[tax$Group == tip, "Species"]
    }
    # get distribution
    tip_distributions[[tip]] = mat_gpl_sum[tip_descendants, "pres" ]
    tip_counts[[tip]] = as.numeric(table(dat [ dat$gene_type == type & dat$taxa %in% tip_descendants, "taxa"]) [ tip_descendants ])
    
  }
  
  # barplots and boxplots
  # presences in extant nodes
  pdf(file = sprintf("%s/phylo_global-type_groups.%s.descendant_dist.pdf", out_fo, type),width = 4, height = 5.5)
  par(mar=c(3.1, 8.1, 4.1, 2.1))
  boxplot(rev(tip_distributions), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", ylim=c(0,100), cex.lab=0.6, cex.axis=0.6)
  stripchart(rev(tip_distributions), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s OG counts in descendant nodes",type), cex.main=0.7)
  
  # counts in extant nodes
  boxplot(rev(tip_counts), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", cex.lab=0.6, cex.axis=0.6, ylim=c(0,350))
  stripchart(rev(tip_counts), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s gene counts in descendant nodes",type), cex.main=0.7)
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-20,20))
  title(sprintf("Gains and losses in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot gain and losses barplots: extant  
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[phyl_data$is_tip,]$gain,
      loss = phyl_data[phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
}






##### plot global cladograms: euk taxonomic groups, per gene type: SUM #####


# prune species tree until groups of interest are tips:
tax_groups = unique(tax$Group)
phyg = phyl
for (taxi in tax_groups) {
  sps_in_tax_group = tax[tax$Group == taxi, "Species"]
  # drop tips one by one, except for first tip
  tip_to_keep = sps_in_tax_group[1]
  tip_to_drop = sps_in_tax_group[-1]
  for (spi in tip_to_drop) {
    phyg = ape::drop.tip(phy = phyg, tip=spi, trim.internal = TRUE)
  }
  # if group contains more than one species, reassign tip name to 
  # its LCA name. Otherwise, keep tip name (e.g. Ttra = Apusozoa)
  if (length(sps_in_tax_group) > 1) {
    phyg$tip.label[phyg$tip.label == tip_to_keep] = taxi
  }
}
# reassing edge lengths (all 1)
phyg$edge.length = rep(1, length(phyg$edge.length))
ape::write.tree(phyg, file = "data/species_tree_pruned.newick")

# dataframe of edges
phyg_edge           = as.data.frame(phyg$edge)
colnames(phyg_edge) = c("edge_start","edge_end")
phyg_edge$ix_edges = as.numeric(rownames(phyg_edge))
phyg_edge$ends_in_tip = phyg_edge$edge_end <= length(phyg$tip.label)
# dataframe of nodes
phyg_nods = data.frame(taxa = c(phyg$tip.label, phyg$node.label))
phyg_nods$edge_end = as.numeric(rownames(phyg_nods))
phyg_nods$is_tip   = phyg_nods$edge_end <= length(phyg$tip.label)
# merge them
phyg_edge = merge(phyg_edge, phyg_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")


# now aggregated cladograms
for (type in gene_types) {
  
  pdf(file = sprintf("%s/phylo_global-sumtype_groups.%s.pdf", out_fo, type),width = 2, height = 3)
  par(mar = c(2, 2, 4.1, 2))
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  # pruned gain los pres matrices
  mai_gain = mat_gain[cla_boo,]
  mai_gain [ mai_gain < 0.1 ] = 0
  mai_loss = mat_loss[cla_boo,]
  mai_loss [ mai_loss < 0.1 ] = 0
  mai_pres = mat_pres[cla_boo,]
  mai_pres [ mai_pres < 0.1 ] = 0
  # aggregate
  mat_gpl_sum$gain = colSums(mai_gain, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mai_loss, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mai_pres, na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE, all.y=FALSE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%.2f")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%.2f")," / -",sprintf(phyl_data$loss, fmt = "%.2f"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %.1f | +%.1f | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %.1f | +%i | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  edgelabels(text=sprintf("+%.1f",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %.1f | +%.1f | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  dev.off()
  
}




##### plot global cladograms: euk taxonomic groups, per gene class #####

# now aggregated cladograms
for (type in gene_classes) {
  
  pdf(file = sprintf("%s/phylo_global-class_groups.%s.pdf", out_fo, type),width = 2, height = 3)
  par(mar = c(2, 2, 4.1, 2))
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_class == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  #   mat_gpl_sum$gain = colSums(mat_gain[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$loss = colSums(mat_loss[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$pres = colSums(mat_pres[cla_boo,], na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE, all.y=FALSE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%#i")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%#i")," / -",sprintf(phyl_data$loss, fmt = "%#i"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  edgelabels(text=sprintf("+%i",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  dev.off()
  
  # plot distribution of OG presence counts in the descendant tips from each "tip" taxonomic group
  tip_distributions = list()
  tip_counts = list()
  for (tip in phyg$tip.label) {
    
    # if display tip is an extant species (will happen for single-species taxon groups):
    if (tip %in% phyl$tip.label) {
      tip_descendants = tip
    } else {
      tip_descendants = tax[tax$Group == tip, "Species"]
    }
    # get distribution
    tip_distributions[[tip]] = mat_gpl_sum[tip_descendants, "pres" ]
    tip_counts[[tip]] = as.numeric(table(dat [ dat$gene_class == type & dat$taxa %in% tip_descendants, "taxa"]) [ tip_descendants ])
  }
  
  # barplots and boxplots
  
  # presences in extant nodes
  pdf(file = sprintf("%s/phylo_global-class_groups.%s.descendant_dist.pdf", out_fo, type),width = 4, height = 5.5)
  par(mar=c(3.1, 8.1, 4.1, 2.1))
  boxplot(rev(tip_distributions), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", ylim=c(0,100), cex.lab=0.6, cex.axis=0.6)
  stripchart(rev(tip_distributions), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s OG counts in descendant nodes",type), cex.main=0.7)
  
  # counts in extant nodes
  boxplot(rev(tip_counts), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", cex.lab=0.6, cex.axis=0.6, ylim=c(0,350))
  stripchart(rev(tip_counts), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s gene counts in descendant nodes",type), cex.main=0.7)
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-20,20))
  title(sprintf("Gains and losses in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot gain and losses barplots: extant  
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[phyl_data$is_tip,]$gain,
      loss = phyl_data[phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
  
}


##### plot global cladograms: euk taxonomic groups, per gene class: SUM #####


# prune species tree until groups of interest are tips:
tax_groups = unique(tax$Group)
phyg = phyl
for (taxi in tax_groups) {
  sps_in_tax_group = tax[tax$Group == taxi, "Species"]
  # drop tips one by one, except for first tip
  tip_to_keep = sps_in_tax_group[1]
  tip_to_drop = sps_in_tax_group[-1]
  for (spi in tip_to_drop) {
    phyg = ape::drop.tip(phy = phyg, tip=spi, trim.internal = TRUE)
  }
  # if group contains more than one species, reassign tip name to 
  # its LCA name. Otherwise, keep tip name (e.g. Ttra = Apusozoa)
  if (length(sps_in_tax_group) > 1) {
    phyg$tip.label[phyg$tip.label == tip_to_keep] = taxi
  }
}
# reassing edge lengths (all 1)
phyg$edge.length = rep(1, length(phyg$edge.length))
ape::write.tree(phyg, file = "data/species_tree_pruned.newick")

# dataframe of edges
phyg_edge           = as.data.frame(phyg$edge)
colnames(phyg_edge) = c("edge_start","edge_end")
phyg_edge$ix_edges = as.numeric(rownames(phyg_edge))
phyg_edge$ends_in_tip = phyg_edge$edge_end <= length(phyg$tip.label)
# dataframe of nodes
phyg_nods = data.frame(taxa = c(phyg$tip.label, phyg$node.label))
phyg_nods$edge_end = as.numeric(rownames(phyg_nods))
phyg_nods$is_tip   = phyg_nods$edge_end <= length(phyg$tip.label)
# merge them
phyg_edge = merge(phyg_edge, phyg_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")


# now aggregated cladograms
for (type in gene_classes) {
  
  pdf(file = sprintf("%s/phylo_global-sumclass_groups.%s.pdf", out_fo, type),width = 2, height = 3)
  par(mar = c(2, 2, 4.1, 2))
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_class == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  # pruned gain los pres matrices
  mai_gain = mat_gain[cla_boo,]
  mai_gain [ mai_gain < 0.1 ] = 0
  mai_loss = mat_loss[cla_boo,]
  mai_loss [ mai_loss < 0.1 ] = 0
  mai_pres = mat_pres[cla_boo,]
  mai_pres [ mai_pres < 0.1 ] = 0
  # aggregate
  mat_gpl_sum$gain = colSums(mai_gain, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mai_loss, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mai_pres, na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE, all.y=FALSE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%.2f")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%.2f")," / -",sprintf(phyl_data$loss, fmt = "%.2f"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %.1f | +%.1f | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %.1f | +%i | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  edgelabels(text=sprintf("+%.1f",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %.1f | +%.1f | -%.1f", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  dev.off()
  
}






##### plot global cladograms: euk taxonomic groups, per gene family #####

# now aggregated cladograms
for (type in unique(mat_summ$og_class)) {
  
  pdf(file = sprintf("%s/phylo_global-family_groups.%s.pdf", out_fo, type),width = 2, height = 3)
  par(mar = c(2, 2, 4.1, 2))
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$og_class == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  #   mat_gpl_sum$gain = colSums(mat_gain[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$loss = colSums(mat_loss[cla_boo,], na.rm = TRUE)
  #   mat_gpl_sum$pres = colSums(mat_pres[cla_boo,], na.rm = TRUE)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = TRUE, all.y=FALSE)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres) / 2)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain) / 2)
  edgelabels(text=paste(phyl_data$taxa,sprintf(phyl_data$pres, fmt = "%#i")), adj = c(0.5,-.2),
             col = alpha("blue",0.4), frame="none", cex=0.4)
  edgelabels(text=paste("+",sprintf(phyl_data$gain, fmt = "%#i")," / -",sprintf(phyl_data$loss, fmt = "%#i"), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.4), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain) / 2)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  edgelabels(text=sprintf("+%i",phyl_data$gain), col = alpha("purple",0.4), frame="none", cex=0.4)
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = "darkgray", root.edge = FALSE, align.tip.label = TRUE, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss) / 2)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss),
    cex.main=0.5, cex.sub=0.5)
  
  # then legend of node size
  node_sizes = c(1,2,5,10,20,50,100)
  plot(y=node_sizes,x=rep(0,length(node_sizes)), cex=sqrt(node_sizes) / 2, log="y", col="springgreen3", pch=16, las=1)
  points(y=5, x=1, cex=sqrt(phyl_data_root$pres) / 2, col=alpha("blue",0.5), pch=16)
  
  dev.off()
  
  # plot distribution of OG presence counts in the descendant tips from each "tip" taxonomic group
  tip_distributions = list()
  tip_counts = list()
  for (tip in phyg$tip.label) {
    
    # if display tip is an extant species (will happen for single-species taxon groups):
    if (tip %in% phyl$tip.label) {
      tip_descendants = tip
    } else {
      tip_descendants = tax[tax$Group == tip, "Species"]
    }
    # get distribution
    tip_distributions[[tip]] = mat_gpl_sum[tip_descendants, "pres" ]
    tip_counts[[tip]] = as.numeric(table(dat [ dat$gene_class == type & dat$taxa %in% tip_descendants, "taxa"]) [ tip_descendants ])
  }
  
  # barplots and boxplots
  
  # presences in extant nodes
  pdf(file = sprintf("%s/phylo_global-family_groups.%s.descendant_dist.pdf", out_fo, type),width = 4, height = 5.5)
  par(mar=c(3.1, 8.1, 4.1, 2.1))
  boxplot(rev(tip_distributions), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", ylim=c(0,100), cex.lab=0.6, cex.axis=0.6)
  stripchart(rev(tip_distributions), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s OG counts in descendant nodes",type), cex.main=0.7)
  
  # counts in extant nodes
  boxplot(rev(tip_counts), horizontal = TRUE, las=1, col = "lightgray", outline = FALSE, border = "gray", cex.lab=0.6, cex.axis=0.6, ylim=c(0,350))
  stripchart(rev(tip_counts), vertical = FALSE, method = "jitter", pch=1, col="blue", add = TRUE, cex = 0.4)
  title(sprintf("%s gene counts in descendant nodes",type), cex.main=0.7)
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-20,20))
  title(sprintf("Gains and losses in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot gain and losses barplots: extant  
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[phyl_data$is_tip,]$gain,
      loss = phyl_data[phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = TRUE, xlab = "Gains or losses", beside = TRUE,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 5.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = TRUE, xlab = "Presence", beside = TRUE,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type),
        cex.main=0.5, cex.sub=0.5)
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
  
}



#### plot individual histories ####

for (type in gene_types) {
  
  print(sprintf("individual histories %s", type)  )
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  
  # list of ogs
  list_ogs = rownames(mat_summ[cla_boo,])
  
  pdf(file = sprintf("%s/histories.%s.pdf", out_fo, type),width = 12, height = 10)
  layout(matrix(1:3, nrow = 1))
  for (ogi in list_ogs) {
    
    ogi_name = mat_summ[ogi,"og_name"]
    ogi_freq = mat_summ[ogi,]$sum_pres_extant
    ogi_fullname = paste(ogi, ogi_name, sep=":")
    ogi_gene_list = dat[dat$orthogroup == ogi_fullname,"gene"]
    
    # find architectures
    ari = arq[arq$gene %in% ogi_gene_list, c("gene","simarq") ]
    if (nrow(ari) == 0) {
      ari = data.frame(frequency = 1, simarq = "absent")
    }
    ari_t = table(ari$simarq)
    
    
    if (ogi_freq > 1) {
      
      ogi_evol = data.frame(row.names = colnames(mat_pres))
      ogi_evol$taxa = rownames(ogi_evol)
      ogi_evol$pres = as.logical(mat_pres[ogi,] > prob_thr)
      ogi_evol$gain = as.logical(mat_gain[ogi,] > prob_thr)
      ogi_evol$loss = as.logical(mat_loss[ogi,] > prob_thr)
      ogi_ancestral_prob = mat_pres[ogi,colnames(mat_pres)[ncol(mat_pres)]]
      
      phyl_ogi_evol = merge(phyl_edge, ogi_evol, by.x = "taxa", by.y = "taxa",all.x = TRUE)
      phyl_ogi_evol = phyl_ogi_evol[order(phyl_ogi_evol$ix_edges),]
      
      phyl_ogi_evol$color = "darkgray"
      phyl_ogi_evol$width = 1
      phyl_ogi_evol$pres [ is.na(phyl_ogi_evol$pres) ] = FALSE
      phyl_ogi_evol[phyl_ogi_evol$pres,"color"] = "lightskyblue"
      phyl_ogi_evol[phyl_ogi_evol$pres,"width"] = 4
      
      # loss nodes
      nodes_w_loss   = as.character(phyl_ogi_evol$taxa[phyl_ogi_evol$loss])
      nodes_w_loss_e = phyl_ogi_evol$ix_edges[phyl_ogi_evol$loss]
      
      # gain node
      nodes_w_gain   = as.character(phyl_ogi_evol$taxa[phyl_ogi_evol$gain])
      nodes_w_gain_e = phyl_ogi_evol$ix_edges[phyl_ogi_evol$gain]
      
      # which tips have ogs?
      tips_w_pres_ix = phyl_ogi_evol[phyl_ogi_evol$pres > 0  & phyl_ogi_evol$edge_end <= length(phyl$tip.label),"edge_end"]
      tips_w_absc_ix = phyl_ogi_evol[phyl_ogi_evol$pres == 0 & phyl_ogi_evol$edge_end <= length(phyl$tip.label),"edge_end"]
      
      # plot phylogeny
      ogi_phylo = plot.phylo(
        phyl, font=1, type="c", label.offset = 4, 
        edge.color = phyl_ogi_evol$color, cex=0.6,
        edge.width = phyl_ogi_evol$width,use.edge.length = FALSE,
        root.edge = TRUE)
      
      # plot presences
      tiplabels(pch = 19, col = "blue",        height = 4, cex = 0.6, adj = +30, tip = tips_w_pres_ix)
      tiplabels(pch = 21, col = "blue", bg=NA, height = 4, cex = 0.6, adj = +30, tip = tips_w_absc_ix)
      
      # plot gain node (if any)
      if (length(nodes_w_gain) > 0) {
        edgelabels(text = paste("+",nodes_w_gain, sep=""), col = "springgreen4", font=2,
                   edge = nodes_w_gain_e, cex=0.6,
                   frame = "none")
      }
      # plot losses (if any)
      if (length(nodes_w_loss_e) > 0) {
        edgelabels(text = paste("!",nodes_w_loss,sep=""), col = "red3", font=2,
                   edge = nodes_w_loss_e, cex=0.6,
                   frame = "none")
      }
      title(main = sprintf(
        "%s\n%s",ogi, ogi_name),
        sub  = sprintf(
          "Gain: %s\nPresence Eukaryota prob = %.2f\nPresent: %i\nLosses: %i", 
          paste(phyl_ogi_evol$taxa[phyl_ogi_evol$gain], collapse = ","), 
          ogi_ancestral_prob,
          sum(phyl_ogi_evol$pres & phyl_ogi_evol$is_tip), 
          sum(phyl_ogi_evol$loss)),
        cex.main=0.7)
      
      # plot architecture heatmap per species
      pie(ari_t, labels = paste(names(ari_t), "\nn=",ari_t, sep=""), cex=0.6)
      
      # plot architecture frequencies
      b=barplot(ari_t, border = NA, xlab = "Frequency",
                col = "lightgray", cex.names = 0.6, horiz = TRUE, width = 1, las=1, ylim = c(0,50))
      text(x=0, y=b, paste(names(ari_t), " | n=",ari_t, sep=""), adj=0, cex=0.6)
      
    }
  }
  dev.off()
  
}






#### main OGs summaries ####

# read taxonomy (again)
taxe = read.table(tax_fn, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
taxe = taxe[match(phyl$tip.label, taxe$Species),]
sps_gap_ixs = c(1, 1 + which(diff(as.numeric(taxe$Group)) != 0)) - 1



for (type in gene_types) {
  
  print(sprintf("presence of main OGs with alikes %s", type)  )
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,] > prob_thr, na.rm = TRUE)
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,] > prob_thr, na.rm = TRUE)
  
  # list of ogs
  list_ogs = rownames(mat_summ[cla_boo,])
  
  pdf(file = sprintf("%s/table_main_OGs.%s.pdf", out_fo, type),width = 15, height = 8)
  layout(matrix(1:3, nrow = 1))
  for (ogi in list_ogs) {
    
    ogi_dom = stringr::str_split(ogi, "\\.", simplify = TRUE)[,1]
    ogi_name = mat_summ[ogi,"og_name"]
    ogi_freq = mat_summ[ogi,]$sum_pres_extant
    ogi_fullname = paste(ogi, ogi_name, sep=":")
    ogi_gene_list = dat[dat$orthogroup == ogi_fullname,"gene"]
    
    # find OGs that are closely related to this one:
    if (ogi_name != "NA" & !grepl("^like:",ogi_name) ) {
      
      ogi_name_items = str_split(ogi_name,"/")[[1]]
      ogi_alike = c()
      for (ogiit in ogi_name_items) {
        ogi_alike = c(ogi_alike,mat_summ[grepl(sprintf("\\b%s\\b",ogiit), mat_summ$og_name) & mat_summ$og_class == ogi_dom & mat_summ$og_id != ogi, "og_id"])
      }
      ogi_alike = unique(ogi_alike)
      # ogi_alike = mat_summ[grepl(ogi_name, mat_summ$og_name) & mat_summ$og_class == ogi_dom & mat_summ$og_id != ogi, "og_id"]
      ogi_alike = c(ogi,ogi_alike)
      ogi_alike_names = unlist(lapply(ogi_alike, function(x) mat_summ[mat_summ$og_id == x ,"og_name_short"]))
      
      ogi_alike_table = mat_pres[ogi_alike,phyl$tip.label]
      rownames(ogi_alike_table) = paste(rownames(ogi_alike_table), ogi_alike_names, sep = " | ")
      ogi_alike_table = rbind(ogi_alike_table, colSums(ogi_alike_table))
      rownames(ogi_alike_table)[nrow(ogi_alike_table)] = "Sum"
      # ogi_alike_table = ogi_alike_table[,order(colSums(ogi_alike_table), decreasing = TRUE)]
      ogi_alike_heatmap = pheatmap(
        ogi_alike_table,  
        color = col_blue(5), breaks = seq(0,5,length.out = 6) - 0.01, 
        cellwidth = 4, cellheight = 4, na_col = "grey", 
        cluster_rows = FALSE, cluster_cols = FALSE, 
        fontsize = 5, 
        gaps_col = sps_gap_ixs, 
        gaps_row = c(1, nrow(ogi_alike_table) - 1),
        legend = FALSE,
        labels_col = taxe$Species.name,
        main = sprintf("%s, %s\n(%i similar OGs)", ogi, ogi_name, length(ogi_alike) - 1),
        border_color = "white")
      
    }    
    
  }
  dev.off()
}



#### ref species summaries ####

# find OGs in each species (for each gene family separately?) and infer their origins

ref_tax_list = c("Hsap","Dmel","Nvec","Aque","Tadh","Metazoa","Cfra","Clim","Cowc","Opisthokonta","Fungi","Scer","Spom","Ncra","Atha","Vcar","Gthe","Ehux","Carmem","Jaklib")
prob_gain_thr = 0.9

mat_pres [ is.na(mat_pres) ] = 0
mat_gain [ is.na(mat_gain) ] = 0

pdf(file = sprintf("%s/reference_taxa_gains.pdf", out_fo, type),width = 6.5, height = 6.5)
parmarori=par()$mar
for (refi in ref_tax_list) {
  
  layout(mat = matrix(1:16, nrow = 4))
  # which nodes are ancestral to the reference node?
  refi_ancestors = c(refi)
  ix_parent = 0
  itertax = refi
  while (!is.null(ix_parent)) {
    ix_parent = phytools::getParent(phyl,node = phyl_nods[phyl_nods$taxa == itertax,"edge_end"])
    nm_parent = as.character(phyl_nods[phyl_nods$edge_end == ix_parent,"taxa"])
    refi_ancestors = c(refi_ancestors, nm_parent)
    itertax = nm_parent
  }
  refi_ancestors = refi_ancestors[refi_ancestors != "Eukaryota"]
  par(mar=parmarori)
  for (type in gene_types) {
    
    # summary
    cla_boo = mat_summ$gene_type == type
    
    # which OGs are present in this taxa?
    ix_present_in_ref = which(mat_pres[cla_boo,refi] >= prob_gain_thr)
    og_present_in_ref = as.character(mat_summ[cla_boo,"og_id"][ix_present_in_ref])
    
    # where were they gained?
    og_gain_in_ref = mat_gain[og_present_in_ref,refi_ancestors]
    og_gain_in_ref_sum = colSums(og_gain_in_ref)
    og_gain_in_ref_sum["Preexistent"] = nrow(og_gain_in_ref) - sum(og_gain_in_ref_sum)
    if (sum(og_gain_in_ref_sum, na.rm = TRUE) > 0) {
      pie(
        og_gain_in_ref_sum,
        col=rainbow(length(og_gain_in_ref_sum), v = 0.9, end = 0.9),
        labels = sprintf("%s\np=%.2f", names(og_gain_in_ref_sum), og_gain_in_ref_sum),
        main = sprintf("Gains of %s %s\nn=%i at p>%.2f", refi, type, nrow(og_gain_in_ref),prob_gain_thr),
        cex.main=0.8,cex=0.6
      )
    } else {
      pie(1)
    }
    
  }
  par(mar=c(0,0,0,0))
  pie(1, col = NA, border = NA, labels=NA)
  legend("topright", legend=names(og_gain_in_ref_sum), fill=rainbow(length(og_gain_in_ref_sum), v = 0.9, end = 0.9), cex=0.4)
  
}
dev.off()





pdf(file = sprintf("%s/reference_taxa_gains_readers_v_all.pdf", out_fo, type),width = 4.5, height = 4.5)
parmarori=par()$mar
layout(mat = matrix(1:12, nrow = 3))
for (refi in ref_tax_list) {
  
  # which nodes are ancestral to the reference node?
  refi_ancestors = c(refi)
  ix_parent = 0
  itertax = refi
  while (!is.null(ix_parent)) {
    ix_parent = phytools::getParent(phyl,node = phyl_nods[phyl_nods$taxa == itertax,"edge_end"])
    nm_parent = as.character(phyl_nods[phyl_nods$edge_end == ix_parent,"taxa"])
    refi_ancestors = c(refi_ancestors, nm_parent)
    itertax = nm_parent
  }
  refi_ancestors = refi_ancestors[refi_ancestors != "Eukaryota"]
  par(mar=parmarori)
  
  
  ## READERS
  # summary
  cla_boo = mat_summ$gene_type == "Readers"
  # which OGs are present in this taxa?
  ix_present_in_ref = which(mat_pres[cla_boo,refi] >= prob_gain_thr)
  og_present_in_ref = as.character(mat_summ[cla_boo,"og_id"][ix_present_in_ref])
  
  # where were they gained?
  og_gain_in_ref = mat_gain[og_present_in_ref,refi_ancestors]
  og_gain_in_ref_sum = colSums(og_gain_in_ref)
  og_gain_in_ref_sum["Preexistent"] = nrow(og_gain_in_ref) - sum(og_gain_in_ref_sum)
  if (sum(og_gain_in_ref_sum, na.rm = TRUE) > 0) {
    pie(
      og_gain_in_ref_sum,
      col=rainbow(length(og_gain_in_ref_sum), v = 0.9, end = 0.9),
      labels = sprintf("%s\np=%.2f", names(og_gain_in_ref_sum), og_gain_in_ref_sum),
      main = sprintf("Gains of %s %s\nn=%i at p>%.2f", refi, "Readers", nrow(og_gain_in_ref),prob_gain_thr),
      cex.main=0.8,cex=0.6
    )
  } else {
    pie(1)
  }
  
  ## CATALYTIC
  # summary
  cla_boo = mat_summ$gene_type %in% c("Acetylase" , "Deacetylase", "Methylase", "Demethylase","Remodeller")
  # which OGs are present in this taxa?
  ix_present_in_ref = which(mat_pres[cla_boo,refi] >= prob_gain_thr)
  og_present_in_ref = as.character(mat_summ[cla_boo,"og_id"][ix_present_in_ref])
  
  # where were they gained?
  og_gain_in_ref = mat_gain[og_present_in_ref,refi_ancestors]
  og_gain_in_ref_sum = colSums(og_gain_in_ref)
  og_gain_in_ref_sum["Preexistent"] = nrow(og_gain_in_ref) - sum(og_gain_in_ref_sum)
  if (sum(og_gain_in_ref_sum, na.rm = TRUE) > 0) {
    pie(
      og_gain_in_ref_sum,
      col=rainbow(length(og_gain_in_ref_sum), v = 0.9, end = 0.9),
      labels = sprintf("%s\np=%.2f", names(og_gain_in_ref_sum), og_gain_in_ref_sum),
      main = sprintf("Gains of %s %s\nn=%i at p>%.2f", refi, "Catalytic", nrow(og_gain_in_ref),prob_gain_thr),
      cex.main=0.8,cex=0.6
    )
  } else {
    pie(1)
  }
  
  par(mar=c(0,0,0,0))
  pie(1, col = NA, border = NA, labels=NA)
  legend("topright", legend=names(og_gain_in_ref_sum), fill=rainbow(length(og_gain_in_ref_sum), v = 0.9, end = 0.9), cex=0.4)
  
}
dev.off()



print("Done!")
