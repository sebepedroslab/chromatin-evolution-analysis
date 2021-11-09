# load libraries
library(ape)
library(scales)
library(phytools)
library(stringr)
library(vioplot)
library(dbscan)


#### Define input ####

domlist=c("SET","BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","HIF-1","Histone",
          "LinkerHistone","PHD","PTIP","SIR2","SNF2_N","TUDOR",
          "zf-CCHH","zf-CXXC","Hist_deacetyl","Acetyltransf_1")



#### plot global cladogram ####

pdf("results_viruses/viral_pairwisetree_radial.pdf", width = 12, height = 12)
for (dom in domlist) {
  
  print(dom)
  
  # list phylogenies for this domain
  phyl_fn = sprintf("results_viruses/pairwise_alignments/vir.%s.diamond.newick", dom)
  
  # load phylogeny
  phyl = read.tree(phyl_fn)
  phyl = midpoint.root(phyl)
  phyl = ladderize(phyl)
  
  # dataframe of edges
  phyl_edge           = as.data.frame(phyl$edge)
  colnames(phyl_edge) = c("edge_start","edge_end")
  phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
  phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
  # dataframe of nodes
  phyl_nods = data.frame(node = c(phyl$tip.label, phyl$node.label))
  phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
  phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)
  
  # taxa in each node
  phyl_nods$taxa = ""
  phyl_nods[phyl_nods$is_tip,"taxa"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "_", simplify = T)[,1]
  # phyl_nods[!phyl_nods$taxa %in% c("arc","bac","vir"),"taxa"] = "euk"
  # gene name
  phyl_nods$gene = ""
  phyl_nods[phyl_nods$is_tip,"gene"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "\\|", simplify = T)[,1]
  
  # merge them
  phyl_edge = merge(phyl_edge, phyl_nods, all.x = T, all.y = F, by.x = "edge_end", by.y = "edge_end",sort = F)
  
  # taxa in each edge
  phyl_edge$taxa = stringr::str_split(phyl_edge$node,pattern = "_", simplify = T)[,1]
  # phyl_edge[!phyl_edge$taxa %in% c("arc","bac","vir"),"taxa"] = "euk"
  
  # add colors
  phyl_edge$color = "slategray3"
  phyl_edge[phyl_edge$taxa %in% c("bac"), "color"] = "blue"
  phyl_edge[phyl_edge$taxa %in% c("arc"), "color"] = "darkturquoise"
  phyl_edge[phyl_edge$taxa %in% c("vir"), "color"] = "deeppink"
  
  # add colors
  phyl_nods$color = "slategray4"
  phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "blue"
  phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = "darkturquoise"
  phyl_nods[phyl_nods$taxa %in% c("vir"), "color"] = "deeppink"
  
  # plot unrooted
  plot.phylo(phyl, font=1, type="u", edge.color = "slategray3", root.edge = T, show.tip.label = F, underscore = T)
  tiplabels(col = phyl_nods$color,frame = "none", pch = 19, cex=0.8)
  tiplabels(text = phyl_nods[phyl_nods$taxa=="vir","gene"],tip = which(phyl_nods$taxa == "vir"), frame="none", col=alpha("deeppink2",0.6), cex=0.4)
  add.scale.bar()
  title(sprintf("%s\n%s\nn=%i sequences",dom,basename(phyl_fn), length(phyl$tip.label)))
  
  legend("topright", col=c("slategray4","blue","darkturquoise","deeppink"), pch=19, legend = c("eukaryotes","bacteria","archaea","viruses"))
}
dev.off()


#### loop distance distribution ####
dist_ta = data.frame()
for (dom in domlist) {
  
  print(dom)
  
  phyl_fn = sprintf("results_viruses/pairwise_alignments/vir.%s.diamond.newick", dom)
  
  # load phylogeny
  phyl = read.tree(phyl_fn)
  
  # indexes
  ix_vir = grep(x = phyl$tip.label, pattern = "^vir_")
  ix_bac = grep(x = phyl$tip.label, pattern = "^bac_")
  ix_arc = grep(x = phyl$tip.label, pattern = "^arc_")
  ix_euk = grep(x = phyl$tip.label, pattern = "^arc_|^vir_|^bac_", invert = T)
  
  # pairwise distances
  dist = ape::cophenetic.phylo(phyl)
  # dist_eps = median(dist)/2
  dist_eps = max(quantile(dist, 0.25),1)
  # dist_eps = 2
  # hist(dist, main=dom, sub=sprintf("eps = %.2f", dist_eps))
  # abline(v=dist_eps, lty=2, col="red")
  
  # find clusters of viral sequences (ignore topology)
  dist_vir = dist[ix_vir, ix_vir]
  dist_d = as.dist(dist_vir)
  dist_k = dbscan(dist_d, eps=dist_eps, minPts = 0)
  
  
  # aggregate distances from viral clusters to other viral clusters
  dist_vir_agg = aggregate(dist_vir, list(dist_k$cluster), median)
  rownames(dist_vir_agg) = dist_vir_agg$Group.1
  dist_vir_agg = dist_vir_agg[,-1]
  dist_vir_agg = t(dist_vir_agg)
  dist_vir_agg = aggregate(dist_vir_agg, list(dist_k$cluster), median)
  rownames(dist_vir_agg) = dist_vir_agg$Group.1
  dist_vir_agg = dist_vir_agg[,-1]
  dist_vir_agg = as.matrix(dist_vir_agg)
  diag(dist_vir_agg) = Inf
  
  # aggregate distances from viral clusters to other cellular sequences
  print(sprintf("%s | eps = %.2f",basename(phyl_fn),dist_eps))
  dist_vir_tocell = dist[ix_vir, -ix_vir]
  dist_vir_tocell_dims = dim(dist_vir_tocell)
  if (!is.null(dist_vir_tocell_dims)) {
    if (dist_vir_tocell_dims[1] > 1 & dist_vir_tocell_dims[2] > 1) {
      # normal scenario
      dist_vic_agg = aggregate(dist_vir_tocell, list(dist_k$cluster), median)
      dist_vic_agg = dist_vic_agg[,-1]
      dist_to_viral = apply(dist_vir_agg, 1, FUN = function(x) min(x))
      dist_to_cells = apply(dist_vic_agg, 1, FUN = function(x) min(x))
      closest_cells = colnames(dist_vic_agg)[apply(dist_vic_agg, 1, FUN = function(x) which.min(x))]
      closest_taxon = stringr::str_split(closest_cells, pattern = "_", simplify = T)[,1]
      closest_taxon[!closest_taxon %in% c("bac","arc")] = "euk"
      # get genes per cluster    
      genes_per_clu = aggregate(rownames(dist_vir), list(dist_k$cluster), paste, collapse=",")[,2]
      numge_per_clu = aggregate(rownames(dist_vir), list(dist_k$cluster), length)[,2]
    } else {
      # in the event of all-virus phylogenies...
      print("Is this an all-virus phylogeny?")
      dist_to_viral = 0
      dist_to_cells = Inf
      closest_taxon = "vir"
      # get genes per cluster    
      genes_per_clu = aggregate(rownames(dist_vir), list(dist_k$cluster), paste, collapse=",")[,2]
      numge_per_clu = aggregate(rownames(dist_vir), list(dist_k$cluster), length)[,2]
    }
  } else {
    # in the event of only one virus within another clade...
    dist_to_viral = Inf
    dist_to_cells = min(dist_vir_tocell)
    closest_cells = names(which.min(dist_vir_tocell))
    closest_taxon = stringr::str_split(closest_cells, pattern = "_", simplify = T)[,1]
    closest_taxon[!closest_taxon %in% c("bac","arc")] = "euk"
    genes_per_clu = rownames(dist)[ix_vir]
    numge_per_clu = length(rownames(dist)[ix_vir])
  }
  
  # relative distance of each viral cluster to other viruses or cells
  dist_t = data.frame(
    dist_to_viral = dist_to_viral,
    dist_to_cells = dist_to_cells,
    closest_taxon = closest_taxon,
    genes_per_clu = genes_per_clu,
    numge_per_clu = numge_per_clu
  )
  dist_t$l2ratio_vir2cell = log2(dist_t$dist_to_viral / dist_t$dist_to_cells)
  dist_t$family = dom
  dist_t$family_phy = basename(phyl_fn)
  
  # store
  dist_ta = rbind(dist_ta, dist_t)
  
  
}


#### plot distance distribution ####
dist_ta$closest_taxon_color = "slategray3"
dist_ta[dist_ta$closest_taxon=="euk","closest_taxon_color"] = "purple"
dist_ta[dist_ta$closest_taxon=="vir","closest_taxon_color"] = "deeppink"
dist_ta[dist_ta$closest_taxon=="bac","closest_taxon_color"] = "blue"
dist_ta[dist_ta$closest_taxon=="arc","closest_taxon_color"] = "darkturquoise"


pdf("results_viruses/viral_distance_distribution.pdf", width = 6, height = 5)
# all families
threshold_def=0.1
order_dist_ta = order(dist_ta$l2ratio_vir2cell)
dist_ta = dist_ta[order_dist_ta,]
plot(dist_ta$l2ratio_vir2cell, col=dist_ta$closest_taxon_color, ylim=c(-3,6), ylab="log2(d to viral / d to cell)")
title(main = "relative distance",
      sub = sprintf("vir-like: %i | cell-like: %i | total: %i",
                    sum(dist_ta$l2ratio_vir2cell < -threshold_def, na.rm = T),
                    sum(dist_ta$l2ratio_vir2cell > threshold_def,  na.rm = T),
                    length(dist_ta$l2ratio_vir2cell)
      )
)
legend("topleft", legend = c("euk","arc","bac"), pch=1, bty = "n", col = c("purple","darkturquoise","blue"))
abline(h=c(0,log2(0.9),log2(1.1)), lty=2, col=c("black","gray","gray"))
# for each family
for (dom in domlist) {
  dist_ti = dist_ta[dist_ta$family == dom,]
  plot(dist_ti$l2ratio_vir2cell, col=dist_ti$closest_taxon_color, ylim=c(-3,6),ylab="log2(d to viral / d to cell)")
  title(main = sprintf("relative distance | %s", dom), 
        sub = sprintf("vir-like: %i | cell-like: %i | total: %i",
                      sum(dist_ti$l2ratio_vir2cell < -threshold_def, na.rm = T),
                      sum(dist_ti$l2ratio_vir2cell > threshold_def,  na.rm = T),
                      length(dist_ti$l2ratio_vir2cell)
        )
  )
  legend("topleft", legend = c("euk","arc","bac"), pch=1, bty = "n", col = c("purple","darkturquoise","blue"))
  abline(h=c(0,-0.1,0.1), lty=2, col=c("black","gray","gray"))
}
dev.off()

