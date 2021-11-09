#### Input ####
# load libraries
library(ape)
library(scales)
library(stringr)
library(pheatmap)
library(ggtree)

gain_fn = "orthogroups_euk.possom.dollo_gain.csv"
loss_fn = "orthogroups_euk.possom.dollo_loss.csv"
pres_fn = "orthogroups_euk.possom.dollo_pres.csv"
summ_fn = "orthogroups_euk.possom.dollo_summary.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
arqs_fn = "architectures.csv"
outp_fn = "results_evol"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"

graphics.off()

#### Define input ####
# load data
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
arq = read.table(arqs_fn, sep = "\t", stringsAsFactors = F, header = F, col.names = c("gene","architecture"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

# fix architecture names
simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
arq$simarq = simple_arqs_vect

gene_counts = table(dat$gene)
genes_repeated = names(gene_counts[gene_counts>1])
ogs_lis_w_dup_genes = dat[dat$gene %in% genes_repeated,"orthogroup"]
dat_dups_full = dat[dat$orthogroup %in% ogs_lis_w_dup_genes,]

# list ogs
ogs_lis = unique(dat_dups_full$orthogroup)
ogs_lis_black = c()

# for (ni in seq_along(ogs_lis)) {
#   
#   ogi = ogs_lis[ni]
#   gli = dat[dat$orthogroup==ogi,"gene"]
#   for (nj in seq_along(ogs_lis)) {
#     
#     if (ni > nj) {
#       
#       ogj = ogs_lis[nj]
#       glj = dat[dat$orthogroup==ogj,"gene"]
#       
#       jac_ij = length(intersect(gli,glj)) / length(union(gli,glj))
#       
#       if (jac_ij > 0.5) {
#         print(sprintf("overlap | jaccard = %.2f | ∩ = %i | u = %i | %s (%i) %s (%i)", 
#                       jac_ij,
#                       length(intersect(gli,glj)),length(union(gli,glj)),
#                       ogi, length(gli), ogj, length(glj)))
#         
#         og_to_remove = ifelse(grepl("Homeobox", ogi), yes = ogj, no= ogi)
#         ogs_lis_black = c(ogs_lis_black, og_to_remove)
#         
#       }
#     }
#   }
# }




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

# define colors for taxa (ancestral and present)
# empty list
taxa_colors = list(gain = c())
# fill extant (all same color)
# for (noi in phyl$tip.label)  { 
#   taxa_colors$gain[noi] = "grey"
# }
# fill ancestral (rainbow)
non = 0
taxa_colors$gain["spsspecific"] = "grey"
for (noi in phyl$node.label) { 
  non=non+1
  taxa_colors$gain[noi] = viridisLite::plasma(length(phyl$node.label),begin = 0.05, end = 0.9)[non]
}


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

# fix OG ids for arqs file
# arq$og_class = stringr::str_split(arq$group,pattern = "\\.", simplify = T)[,1]
# arq$og_id    = stringr::str_split(arq$group,pattern = ":", simplify = T)[,1]


# load gains, losses and ancestral presences
mat_gain = read.table(gain_fn, sep="\t", header = T, row.names = 1)
mat_loss = read.table(loss_fn, sep="\t", header = T, row.names = 1)
mat_pres = read.table(pres_fn, sep="\t", header = T, row.names = 1)

# add num species presence (binarised)
mat_summ$n_presence_bin = apply(mat_pres[,phyl$tip.label], MARGIN = 1, function(x) sum(x>0))

# simplify OG names
rownames(mat_gain) = mat_summ$og_id
rownames(mat_loss) = mat_summ$og_id
rownames(mat_pres) = mat_summ$og_id

# blacklist OGs:
# remove singletons (only present in one species)? -> no
# ogs_lis_black = c(ogs_lis_black, rownames(mat_summ[mat_summ$n_presence_bin==1,]))
# remove OGs where number of parallel losses exceeds number of actual presences
# ogs_lis_black = c(ogs_lis_black, rownames(mat_summ[mat_summ$n_presence_bin == 2 & mat_summ$n_losses > mat_summ$n_presence_bin*5,])   )

# remove blacklisted OGs
print(sprintf("Remove %i OGs", length(ogs_lis_black)))
print(ogs_lis_black)
mat_summ = mat_summ[!(rownames(mat_summ) %in% ogs_lis_black),]
mat_gain = mat_gain[!(rownames(mat_gain) %in% ogs_lis_black),]
mat_loss = mat_loss[!(rownames(mat_loss) %in% ogs_lis_black),]
mat_pres = mat_pres[!(rownames(mat_pres) %in% ogs_lis_black),]

# order taxa according to phylogeny?
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]
mat_loss = mat_loss[,taxa_lis]

# create factors for taxa and families
mat_summ$gain = factor(mat_summ$gain, levels = taxa_lis)
mat_summ$og_class = factor(mat_summ$og_class, levels = unique(mat_summ$og_class))

# add classes to per-gene table
dat = merge(dat, mat_summ, by.x="orthogroup", by.y="orthogroup", all.x = T, all.y = F)
# add species 
dat$taxa = stringr::str_split(dat$gene,pattern = "_", simplify = T)[,1]
dat$taxa = factor(dat$taxa, levels = taxa_lis)
# remove blacklisted ogs
dat = dat[!(dat$og_class %in% ogs_lis_black), ]


# define colors heatmap
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","dodgerblue4"))

# list of gene types
gene_types = unique(as.character(mat_summ$gene_type))
gene_classes = unique(as.character(mat_summ$gene_class))




##### plot global cladograms: all species #####


# first, plot model phylogenies
pdf(file = sprintf("%s/phylo_global_sps.00.pdf", outp_fn),width = 4, height = 8)
ggtree(phyl, ladderize = F, layout = "rectangular",branch.length = "none", lwd=0.5, color="gray") +
  geom_nodelab(color="red", size=1.5) +
  geom_tiplab(offset=0, size=1.5) +
  ggplot2::scale_y_reverse()
dev.off()

# now aggregated cladograms
for (type in gene_types) {
  
  pdf(file = sprintf("%s/phylo_global-type_sps.%s.pdf", outp_fn, type),width = 8, height = 20)
  print(sprintf("Plot aggregated counts %s | species", type))
  
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,])
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,])
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,]>0)
  
  # add info to plot
  phyl_data = merge(phyl_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = T)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres)/4)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3)," / -",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then gains
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3),sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then losses
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("-",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
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
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
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
    phyg = ape::drop.tip(phy = phyg, tip=spi, trim.internal = T)
  }
  # if group contains more than one species, reassign tip name to 
  # its LCA name. Otherwise, keep tip name (e.g. Ttra = Apusozoa)
  if (length(sps_in_tax_group) > 1) {
    phyg$tip.label[phyg$tip.label==tip_to_keep] = taxi
  }
}
# reassing edge lengths (all 1)
phyg$edge.length = rep(1, length(phyg$edge.length))

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
phyg_edge = merge(phyg_edge, phyg_nods, all.x = T, all.y = T, by.x = "edge_end", by.y = "edge_end")


# first, plot model phylogenies
pdf(file = sprintf("%s/phylo_global_groups.00.pdf", outp_fn),width = 3.5, height = 4.5)
ggtree(phyg, ladderize = F, layout = "rectangular",branch.length = "none", lwd=0.5, color="gray") +
  geom_nodelab(color="red", size=1.5) +
  geom_tiplab(offset=0, size=1.5) +
  ggplot2::scale_y_reverse()
dev.off()

# now aggregated cladograms
for (type in gene_types) {
  
  pdf(file = sprintf("%s/phylo_global-type_groups.%s.pdf", outp_fn, type),width = 4, height = 6)
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,])
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,])
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,]>0)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = T, all.y=F)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres)/4)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3)," / -",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3),sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("-",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
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
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
  # plot distribution of OG presence counts in the descendant tips from each "tip" taxonomic group
  tip_distributions = list()
  for (tip in phyg$tip.label) {
    
    # if display tip is an extant species (will happen for single-species taxon groups):
    if (tip %in% phyl$tip.label) {
      tip_descendants = tip
    } else {
      tip_descendants = tax[tax$Group == tip, "Species"]
    }
    # get distribution
    tip_distributions[[tip]] = mat_gpl_sum[tip_descendants, "pres" ]
    
  }
  pdf(file = sprintf("%s/phylo_global-type_groups.%s.descendant_dist.pdf", outp_fn, type),width = 4, height = 8)
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  boxplot(rev(tip_distributions), horizontal = T, las=1, col = "lightgray", outline = F, border = "gray", ylim=c(0,100))
  stripchart(rev(tip_distributions), vertical = F, method = "jitter", pch=1, col="blue", add = T, cex = 0.8)
  title(sprintf("Presence in descendant nodes, %s",type))
  dev.off()
  
}








##### plot global cladograms: euk taxonomic groups, per gene class #####

# now aggregated cladograms
for (type in gene_classes) {
  
  pdf(file = sprintf("%s/phylo_global-class_groups.%s.pdf", outp_fn, type),width = 4, height = 6)
  print(sprintf("Plot aggregated counts %s | groups", type))
  
  cla_boo = mat_summ$gene_class == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,])
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,])
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,]>0)
  
  # add info to plot
  phyl_data = merge(phyg_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = T, all.y=F)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # store root node in another table
  phyl_data_root = phyl_data[is.na(phyl_data$ix_edges),]
  phyl_data = phyl_data[!is.na(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("cyan2",0.2), cex=sqrt(phyl_data$pres)/4)
  edgelabels(pch=16, col = alpha("springgreen3",0.8), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3)," / -",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Presence and gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then gains
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("springgreen3",0.5), cex=sqrt(phyl_data$gain)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("+",signif(phyl_data$gain,3),sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Gains %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # then losses
  plot.phylo(phyg, font=1, type="phylogram", label.offset = 0, edge.color = phyl_data$color, root.edge = F, align.tip.label = T, cex=0.4)
  edgelabels(pch=16, col = alpha("darkorange",0.5), cex=sqrt(phyl_data$loss)/4)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.4)
  edgelabels(text=paste("-",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.4)
  title(
    main=sprintf("Losses %s",type),
    sub=sprintf("Root: %s | %i | +%i | -%i", phyl_data_root$taxa, phyl_data_root$pres, phyl_data_root$gain, phyl_data_root$loss))
  
  # plot gain and losses barplots: ancestral
  gldat = as.matrix(
    t(data.frame(
      gain=phyl_data[!phyl_data$is_tip,]$gain,
      loss = phyl_data[!phyl_data$is_tip,]$loss * -1, 
      row.names = phyl_data[!phyl_data$is_tip,]$taxa)))
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(gldat,cex.names = 0.5,col = c("springgreen3","darkorange"),
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
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
          las=1, horiz = T, xlab = "Gains or losses", beside = T,border = NA,
          xlim = c(-50,50))
  title(sprintf("Gains and losses in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: extant  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[phyl_data$is_tip,]$taxa,"n =",phyl_data[phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in extant nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  # plot presence barplots: ancestral  
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  barplot(phyl_data[!phyl_data$is_tip,]$pres,
          names.arg = paste(phyl_data[!phyl_data$is_tip,]$taxa,"n =",phyl_data[!phyl_data$is_tip,]$pres),
          cex.names = 0.5,col = c("deepskyblue"),
          las=1, horiz = T, xlab = "Presence", beside = T,border = NA,
          xlim = c(0,100))
  title(sprintf("Presence in ancestral nodes, %s",type))
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
  
  # plot distribution of OG presence counts in the descendant tips from each "tip" taxonomic group
  tip_distributions = list()
  for (tip in phyg$tip.label) {
    
    # if display tip is an extant species (will happen for single-species taxon groups):
    if (tip %in% phyl$tip.label) {
      tip_descendants = tip
    } else {
      tip_descendants = tax[tax$Group == tip, "Species"]
    }
    # get distribution
    tip_distributions[[tip]] = mat_gpl_sum[tip_descendants, "pres" ]
    
  }
  pdf(file = sprintf("%s/phylo_global-class_groups.%s.descendant_dist.pdf", outp_fn, type),width = 4, height = 8)
  par(mar=c(3.1, 10.1, 4.1, 2.1))
  boxplot(rev(tip_distributions), horizontal = T, las=1, col = "lightgray", outline = F, border = "gray", ylim=c(0,100))
  stripchart(rev(tip_distributions), vertical = F, method = "jitter", pch=1, col="blue", add = T, cex = 0.8)
  title(sprintf("Presence in descendant nodes, %s",type))
  dev.off()
  
  
}






#### plot individual histories ####

for (type in gene_types) {
  
  print(sprintf("individual histories %s", type)  )
  cla_boo = mat_summ$gene_type == type
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,])
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,])
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,]>0)
  
  # list of ogs
  list_ogs = rownames(mat_summ[cla_boo,])
  
  pdf(file = sprintf("%s/phylo_histories.%s.pdf", outp_fn, type),width = 12, height = 10)
  layout(matrix(1:3, nrow = 1))
  for (ogi in list_ogs) {
    
    ogi_name = mat_summ[ogi,"og_name"]
    ogi_freq = mat_summ[ogi,]$n_presences
    ogi_fullname = paste(ogi, ogi_name, sep=":")
    ogi_gene_list = dat[dat$orthogroup == ogi_fullname,"gene"]
    
    # find architectures
    ari = arq[arq$gene %in% ogi_gene_list, c("gene","simarq") ]
    if (nrow(ari) == 0) {
      ari = data.frame(frequency = 1, simarq = "absent")
    }
    ari_t = table(ari$simarq)
    
    
    if (ogi_freq>1) {
      
      ogi_evol = data.frame(row.names = colnames(mat_pres))
      ogi_evol$taxa = rownames(ogi_evol)
      ogi_evol$pres = as.logical(mat_pres[ogi,] > 0)
      ogi_evol$gain = as.logical(mat_gain[ogi,] > 0)
      ogi_evol$loss = as.logical(mat_loss[ogi,] > 0)
      
      phyl_ogi_evol = merge(phyl_edge, ogi_evol, by.x = "taxa", by.y = "taxa",all.x = T)
      phyl_ogi_evol = phyl_ogi_evol[order(phyl_ogi_evol$ix_edges),]
      
      phyl_ogi_evol$color = "slategray"
      phyl_ogi_evol$width = 1
      phyl_ogi_evol[phyl_ogi_evol$pres,"color"] = "lightskyblue"
      phyl_ogi_evol[phyl_ogi_evol$pres,"width"] = 4
      
      # loss nodes
      nodes_w_loss   = as.character(phyl_ogi_evol$taxa[phyl_ogi_evol$loss])
      nodes_w_loss_e = phyl_ogi_evol$ix_edges[phyl_ogi_evol$loss]
      
      # gain node
      nodes_w_gain   = as.character(phyl_ogi_evol$taxa[phyl_ogi_evol$gain])
      nodes_w_gain_e = phyl_ogi_evol$ix_edges[phyl_ogi_evol$gain]
      
      # which tips have ogs?
      tips_w_pres_ix = phyl_ogi_evol[phyl_ogi_evol$pres>0  & phyl_ogi_evol$edge_end<=length(phyl$tip.label),"edge_end"]
      tips_w_absc_ix = phyl_ogi_evol[phyl_ogi_evol$pres==0 & phyl_ogi_evol$edge_end<=length(phyl$tip.label),"edge_end"]
      
      plot.phylo(phyl, font=1, type="cladogram", label.offset = 4, 
                 edge.color = phyl_ogi_evol$color, cex=0.6,
                 edge.width = phyl_ogi_evol$width,use.edge.length = F,
                 root.edge = T)
      # plot presences
      tiplabels(pch = 19, col = "blue",        height = 4, cex = 0.9, adj = +30, tip = tips_w_pres_ix)
      tiplabels(pch = 21, col = "blue", bg=NA, height = 4, cex = 0.9, adj = +30, tip = tips_w_absc_ix)
      # plot gain node
      edgelabels(text = paste("+",nodes_w_gain, sep=""), col = "springgreen4", font=2,
                 edge = nodes_w_gain_e, cex=0.6,
                 frame = "none")
      # plot losses (if any)
      if (length(nodes_w_loss_e) > 0) {
        edgelabels(text = paste("!",nodes_w_loss,sep=""), col = "red3", font=2,
                   edge = nodes_w_loss_e, cex=0.6,
                   frame = "none")
      }
      title(main = sprintf("%s\n%s",ogi, ogi_name),
            sub  = sprintf("Gain: %s\nPresent: %i\nLosses: %i", 
                           phyl_ogi_evol$taxa[phyl_ogi_evol$gain], 
                           sum(phyl_ogi_evol$pres & phyl_ogi_evol$is_tip), 
                           sum(phyl_ogi_evol$loss)),
            cex.main=0.7)
      
      # plot architecture heatmap per species
      pie(ari_t, labels = paste(names(ari_t), "\nn=",ari_t, sep=""), cex=0.6)
      
      # plot architecture frequencies
      b=barplot(ari_t, border = NA, xlab = "Frequency",
                col = "lightgray", cex.names = 0.6, horiz = T, width = 1, las=1, ylim = c(0,50))
      text(x=0, y=b, paste(names(ari_t), " | n=",ari_t, sep=""), adj=0, cex=0.8)
      
    }
  }
  dev.off()
  
}


#### OG presence heatmaps ####

for (type in gene_types) {
  
  pdf(file = sprintf("%s/heatmap_presence.%s.pdf", outp_fn, type),width = 12, height = 36)
  
  print(sprintf("heatmaps %s", type)  )
  
  list_fams_in_type = unique(as.character(mat_summ[mat_summ$gene_type == type,"og_class"]))
  
  for (fam in list_fams_in_type) {
    
    cla_boo = mat_summ$og_class == fam
    mat_pres_cla = mat_pres[cla_boo,]
    mat_pres_ids = mat_summ$og_id_short[cla_boo]
    
    # order according to number of presences
    mat_pres_ord = order(rowSums(mat_pres_cla>0),decreasing = T)
    mat_pres_cla = mat_pres_cla[mat_pres_ord,]
    mat_pres_ids = mat_pres_ids[mat_pres_ord]
    
    mat_pres_cla_split = split(mat_pres_cla, rep(1:ceiling(nrow(mat_pres_cla)/100), each=100, length.out=nrow(mat_pres_cla)))
    mat_pres_ids_split = split(mat_pres_ids, rep(1:ceiling(nrow(mat_pres_cla)/100), each=100, length.out=nrow(mat_pres_cla)))
    
    mat_gain_cols = mat_summ[rownames(mat_pres_cla),"gain", drop=F]
    mat_gain_cols$gain = as.character(mat_gain_cols$gain)
    mat_gain_cols[mat_gain_cols$gain %in% phyl$tip.label,"gain"] = "spsspecific"
    
    non=0
    for (mat_pres_cla_i in mat_pres_cla_split) {
      
      non=non+1
      pheatmap(t(mat_pres_cla_i), color = col_blue(20), breaks = seq(0,1,length.out = 21), gaps_row = length(phyl$tip.label),
               cellwidth = 6, cellheight = 6, na_col = "grey", cluster_rows = F, cluster_cols = F,
               fontsize = 6,main = sprintf("%s presence (%i/%i)", fam, non, length(mat_pres_cla_split)),
               annotation_col = mat_gain_cols,
               annotation_colors = taxa_colors,
               labels_col = mat_pres_ids_split[[non]],
               border_color = "white", display_numbers = F)
      
    }
    
  }
  
  dev.off()
  
  
}


print("All done!")