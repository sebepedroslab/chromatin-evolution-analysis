# load libraries
library(ape)
library(scales)
library(stringr)
library(pheatmap)

outp_fn = "results_evol/euk"
gain_fn = "orthogroups_euk.dollo_gain.csv"
loss_fn = "orthogroups_euk.dollo_loss.csv"
pres_fn = "orthogroups_euk.dollo_pres.csv"
summ_fn = "orthogroups_euk.dollo_summary.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "species_tree.newick"

graphics.off()

#### Define input ####

# load data
dat = read.table(dat_fn, header = T, stringsAsFactors = F)

gene_counts = table(dat$gene)
genes_repeated = names(gene_counts[gene_counts>1])
ogs_lis_w_dup_genes = dat[dat$gene %in% genes_repeated,"orthogroup"]
dat_dups_full = dat[dat$orthogroup %in% ogs_lis_w_dup_genes,]

# list ogs
ogs_lis = unique(dat_dups_full$orthogroup)
ogs_lis_black = c()

for (ni in seq_along(ogs_lis)) {
  
  ogi = ogs_lis[ni]
  gli = dat[dat$orthogroup==ogi,"gene"]
  for (nj in seq_along(ogs_lis)) {
    
    if (ni > nj) {
      
      ogj = ogs_lis[nj]
      glj = dat[dat$orthogroup==ogj,"gene"]
      
      jac_ij = length(intersect(gli,glj)) / length(union(gli,glj))
      
      if (jac_ij > 0.5) {
        print(sprintf("overlap | jaccard = %.2f | ∩ = %i | u = %i | %s (%i) %s (%i)", 
                      jac_ij,
                      length(intersect(gli,glj)),length(union(gli,glj)),
                      ogi, length(gli), ogj, length(glj)))
        
        og_to_remove = ifelse(grepl("Homeobox", ogi), yes = ogj, no= ogi)
        ogs_lis_black = c(ogs_lis_black, og_to_remove)
        
      }
    }
  }
}

# input files
phyl = read.tree(file = phyl_fn)

# summary
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


# load gains, losses and ancestral presences
mat_gain = read.table(gain_fn, sep="\t", header = T, row.names = 1)
mat_loss = read.table(loss_fn, sep="\t", header = T, row.names = 1)
mat_pres = read.table(pres_fn, sep="\t", header = T, row.names = 1)

# remove blacklisted OGs (redundant with others)
mat_gain = mat_gain[!(rownames(mat_gain) %in% ogs_lis_black),]
mat_loss = mat_loss[!(rownames(mat_loss) %in% ogs_lis_black),]
mat_pres = mat_pres[!(rownames(mat_pres) %in% ogs_lis_black),]

# simplify OG names
rownames(mat_gain) = mat_summ$og_id
rownames(mat_loss) = mat_summ$og_id
rownames(mat_pres) = mat_summ$og_id
rownames(mat_summ) = mat_summ$og_id


# order taxa according to phylogeny?
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]
mat_loss = mat_loss[,taxa_lis]

# create factors for taxa and classes
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
col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","dodgerblue4"))

# define colors for taxa (ancestral and present)
# empty list
taxa_colors = list(gain = c())
# fill extant (all same color)
for (noi in phyl$tip.label)  { 
  taxa_colors$gain[noi] = "grey"
}
# fill ancestral (rainbow)
non = 0
for (noi in phyl$node.label) { 
  non=non+1
  taxa_colors$gain[noi] = viridisLite::viridis(length(phyl$node.label),begin = 0.05, end = 0.9)[non]
}

#### plot global cladogram ####

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

# summarise gains, losses and presences per node
mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
mat_gpl_sum$taxa = rownames(mat_gpl_sum)
mat_gpl_sum$gain = colSums(mat_gain)
mat_gpl_sum$loss = colSums(mat_loss)
mat_gpl_sum$pres = colSums(mat_pres>0)

# add info to plot
phyl_data = merge(phyl_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = T)
phyl_data = phyl_data[order(phyl_data$ix_edges),]

# add colors
phyl_data$color = "slategray"
phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"

# plot phylograms
pdf(file = sprintf("%s.phylo_phylogram.pdf", outp_fn),width = 9, height = 18)

# first, global
plot.phylo(phyl, font=1, type="phylogram", label.offset = 10, edge.color = phyl_data$color, root.edge = T)
edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
           col = alpha("blue",0.8), frame="none", cex=0.7)
edgelabels(text=paste("+",signif(phyl_data$gain,3)," / -",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
           col = alpha("purple",0.8), frame="none", cex=0.7)
title(sprintf("Gains and losses %s","general"))

# second, per-class
# list of classes
list_class = as.character(unique(mat_summ$og_class))
for (cla in list_class) {
  
  cla_boo = mat_summ$og_class == cla
  
  # summarise gains, losses and presences per node
  mat_gpl_sum = data.frame(row.names = colnames(mat_gain))
  mat_gpl_sum$taxa = rownames(mat_gpl_sum)
  mat_gpl_sum$gain = colSums(mat_gain[cla_boo,])
  mat_gpl_sum$loss = colSums(mat_loss[cla_boo,])
  mat_gpl_sum$pres = colSums(mat_pres[cla_boo,]>0)
  
  # add info to plot
  phyl_data = merge(phyl_edge, mat_gpl_sum, by.x = "taxa", by.y = "taxa",all.x = T)
  phyl_data = phyl_data[order(phyl_data$ix_edges),]
  
  # add colors
  phyl_data$color = "slategray"
  phyl_data[phyl_data$gain > phyl_data$loss,"color"] = "springgreen4"
  phyl_data[phyl_data$gain < phyl_data$loss,"color"] = "deeppink3"
  
  # first, global
  plot.phylo(phyl, font=1, type="phylogram", label.offset = 10, edge.color = phyl_data$color, root.edge = T)
  edgelabels(text=paste(phyl_data$taxa,signif(phyl_data$pres,3)), adj = c(0.5,-.2),
             col = alpha("blue",0.8), frame="none", cex=0.7)
  edgelabels(text=paste("+",signif(phyl_data$gain,3)," / -",signif(phyl_data$loss,3), sep=""), adj = c(0.5,1.2),
             col = alpha("purple",0.8), frame="none", cex=0.7)
  title(sprintf("Gains and losses %s",cla))
  
  
  
}
dev.off()


#### plot cladograms per OG ####

# list of ogs
list_ogs = rownames(mat_summ)

pdf(file = sprintf("%s.phylo_individual_histories.pdf", outp_fn),width = 5, height = 15)
for (ogi in list_ogs) {
  
  ogi_name = mat_summ[ogi,"og_name"]
  ogi_freq = mat_summ[ogi,]$n_presences
  
  if (ogi_freq>1) { 
    
    print(paste(ogi, ogi_name, ogi_freq))
    
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
    
    plot.phylo(phyl, font=1, type="phylogram", label.offset = 4, 
               edge.color = phyl_ogi_evol$color, cex=0.7,
               edge.width = phyl_ogi_evol$width,
               root.edge = T)
    # plot presences
    tiplabels(pch = 19, col = "blue",        height = 4, cex = 1, adj = +30, tip = tips_w_pres_ix)
    tiplabels(pch = 21, col = "blue", bg=NA, height = 4, cex = 1, adj = +30, tip = tips_w_absc_ix)
    # plot gain node
    edgelabels(text = paste("+",nodes_w_gain, sep=""), col = "springgreen4", font=2,
               edge = nodes_w_gain_e, cex=0.7,
               frame = "none")
    # plot losses (if any)
    if (length(nodes_w_loss_e) > 0) {
      edgelabels(text = paste("!",nodes_w_loss,sep=""), col = "red3", font=2,
                 edge = nodes_w_loss_e, cex=0.7,
                 frame = "none")
    }
    title(main = sprintf("%s\n%s",ogi, ogi_name),
          sub  = sprintf("Gain: %s\nPresent: %i\nLosses: %i", 
                         phyl_ogi_evol$taxa[phyl_ogi_evol$gain], 
                         sum(phyl_ogi_evol$pres & phyl_ogi_evol$is_tip), 
                         sum(phyl_ogi_evol$loss)),
          cex.main=0.7)
  }
}
dev.off()



#### OG presence heatmaps ####

pdf(file = sprintf("%s.heatmap_presence.pdf", outp_fn),width = 12, height = 36)
for (cla in list_class) {
  
  print(cla)
  cla_boo = mat_summ$og_class == cla
  mat_pres_cla = mat_pres[cla_boo,]
  mat_pres_ids = mat_summ$og_id_short[cla_boo]
  
  # order according to number of presences
  mat_pres_ord = order(rowSums(mat_pres_cla>0),decreasing = T)
  mat_pres_cla = mat_pres_cla[mat_pres_ord,]
  mat_pres_ids = mat_pres_ids[mat_pres_ord]
  
  mat_pres_cla_split = split(mat_pres_cla, rep(1:ceiling(nrow(mat_pres_cla)/100), each=100, length.out=nrow(mat_pres_cla)))
  mat_pres_ids_split = split(mat_pres_ids, rep(1:ceiling(nrow(mat_pres_cla)/100), each=100, length.out=nrow(mat_pres_cla)))
  
  non=0
  for (mat_pres_cla_i in mat_pres_cla_split) {
    
    non=non+1
    pheatmap(t(mat_pres_cla_i), color = col_blue(20), breaks = seq(0,1,length.out = 21), gaps_row = length(phyl$tip.label),
             cellwidth = 6, cellheight = 6, na_col = "grey", cluster_rows = F, cluster_cols = F,
             fontsize = 6,main = sprintf("%s presence (%i/%i)", cla, non, length(mat_pres_cla_split)),
             annotation_col = mat_summ[,"gain", drop=F],
             annotation_colors = taxa_colors,
             labels_col = mat_pres_ids_split[[non]],
             border_color = "white", display_numbers = F)
    
  }
  
  
}
dev.off()



#### OG gains, presence, loss, etc. per node and type barplots ####

# gains
pdf(file = sprintf("%s.pie_gains_per_type_per_node.pdf", outp_fn),width = 4, height = 4)
for (tax in taxa_lis) {
  
  tax_boo_gain = mat_gain[,tax] > 0
  mat_summ_tax = mat_summ[tax_boo_gain,]
  mat_summ_tax_tab = table(mat_summ_tax$og_class)
  mat_summ_tax_mat = matrix(mat_summ_tax_tab)
  rownames(mat_summ_tax_mat) = names(mat_summ_tax_tab)
  
  if (sum(mat_summ_tax_tab)>0) {
    pie(mat_summ_tax_tab, 
        labels = paste(names(mat_summ_tax_tab)," N=", mat_summ_tax_tab, sep=""),
        col=rainbow(nlevels(mat_summ_tax$og_class), v = 0.8),
        main=sprintf("gains %s N=%i",tax,sum(mat_summ_tax_tab)))
  }
}
dev.off()

# presence
pdf(file = sprintf("%s.barplot_pres_per_type_per_node.pdf", outp_fn),width = 4, height = 2)
for (tax in taxa_lis) {
  
  tax_boo_gain = mat_pres[,tax] > 0
  mat_summ_tax = mat_summ[tax_boo_gain,]
  mat_summ_tax_tab = table(mat_summ_tax$og_class)
  mat_summ_tax_mat = matrix(mat_summ_tax_tab)
  rownames(mat_summ_tax_mat) = names(mat_summ_tax_tab)
  
  barplot(mat_summ_tax_mat, horiz = T, 
          col=rainbow(50, v = 0.8),
          main=tax)
}
dev.off()



# presence of genes and OGs per species
gene_og_tab = table(dat$taxa)

# num genes
gene_og_mat = matrix(gene_og_tab)
rownames(gene_og_mat) = names(gene_og_tab)
colnames(gene_og_mat)[1] = "genes"

# add num duplications within each og
mat_num_duplications = mat_pres-1
mat_num_duplications[mat_num_duplications<0] = 0
gene_og_mat = cbind(gene_og_mat, colSums(mat_num_duplications))
colnames(gene_og_mat)[2] = "intra-OG paralogs"

# add num OGs
gene_og_mat = cbind(gene_og_mat, colSums(mat_pres>0))
colnames(gene_og_mat)[3] = "OGs"

# add num OGs with duplication
gene_og_mat = cbind(gene_og_mat, colSums(mat_pres>1))
colnames(gene_og_mat)[4] = "OGs with internal paralogs"

# transpose
gene_og_mat = t(gene_og_mat)

pdf(file = sprintf("%s.barplot_pres_per_node.pdf", outp_fn),width = 4, height = 20)
b=barplot(gene_og_mat,horiz = T, col = c("violet","purple","lightblue","blue"), las=1, beside = T, 
          xlab = "Counts",
          cex.names = 0.6, 
          legend=rownames(gene_og_mat), args.legend = list(cex=0.6))
text(x=max(gene_og_mat),b,
     labels = gene_og_mat, 
     col=c("orange","red"),
     pos=2,cex=0.4)
dev.off()


# legend with per-type colors
pdf(file = sprintf("%s.barplot_legend_per_type_per_node.pdf", outp_fn),width = 4, height = 20)
mat_summ_tax_mat[,1] = 0
barplot(mat_summ_tax_mat, horiz = T, 
        col=rainbow(50, v = 0.8),border = NA,
        main=tax, legend=rownames(mat_summ_tax_mat))
dev.off()


#### OG gain nodes per species ####

pdf(file = sprintf("%s.heatmap_gains_per_node.pdf", outp_fn),width = 4, height = 8)
for (tax in taxa_lis) {
  
  tax_boo_pres = mat_pres[,tax] > 0
  mat_summ_tax = mat_summ[tax_boo_pres,]
  
  # gains per node among ogs present in this taxon
  mat_summ_tax_gain = as.matrix(ftable(xtabs(formula = ~ og_class+gain , data=mat_summ_tax)))
  # present in this taxon
  mat_summ_tax_pres = as.vector(xtabs(formula = ~ og_class , data=mat_summ_tax))
  mat_summ_tax_gain_frac = mat_summ_tax_gain/mat_summ_tax_pres
  
  # remove absent nodes
  mat_summ_tax_gain_frac = mat_summ_tax_gain_frac [, colSums(mat_summ_tax_gain) > 0 ]
  
  pheatmap(mat_summ_tax_gain_frac, color = col_blue(20), breaks = seq(0,1,length.out = 21), 
           cellwidth = 10, cellheight = 6, na_col = "grey", cluster_rows = F, cluster_cols = F,
           number_format = "%.2f", number_color = "red",fontsize_number = 4,
           fontsize = 6,main = sprintf("TFs in %s: OGs gained per node", tax),
           border_color = "white", display_numbers = T)
  
}
dev.off()


stop("arra")
