# load libraries
library(viridis)
library(ape)
library(scales)
library(phytools)
library(adephylo)
library(stringr)
library(data.table)
library(pheatmap)


#### Define input ####

# domains with viral hits
domlist=c("Acetyltransf_1","GNAT_acetyltr_2","MOZ_SAS","Hist_deacetyl","SIR2","DOT1","SET","SNF2_N","ASF1_hist_chap")
# domlist=c("Hist_deacetyl","GNAT_acetyltr_2","DOT1")

# tables
vtax_fn = "../data/seq_Archaea.taxa.csv.gz"    # taxon for each viral seq
taxne_fn = "../data/taxonomy_ranked_wg.tsv.gz" # taxonomy for non-euks
ort_fn = "orthogroups_euk.csv"
cla_fn = "../data/gene_families_hmm.csv"
arq_fn = "gene_sequences/noneuk_arc.pfamscan_archs.csv"
tre_fo = "results_preeuk/trees_ogs/"

# read data
# taxonomy of viral sequences
vtax = fread(cmd=sprintf("zcat %s", vtax_fn), data.table = F, header = F, stringsAsFactors = F)
colnames(vtax) = c("gene", "taxon")
# taxonomy of species
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=F, stringsAsFactors=F)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")
taxne$phylum = gsub("Candidatus ", "",taxne$phylum)
# subset to viral sequences for speed
# taxne = taxne[taxne$superkingdom == "Viruses",]
# orthology, gene fam classification, architectures
ort = read.table(ort_fn, sep="\t", header = T, stringsAsFactors = F)
cla = read.table(
  cla_fn, header = F, sep = "\t", 
  col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
  stringsAsFactors = F)
arq = read.table(arq_fn, header = F, col.names = c("gene","architecture"), sep="\t")

# simplify architectures (collapse **consecutive** repeats)
simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
arq$simarq = simple_arqs_vect

# ordered list of domains
domlist = domlist[match(cla$gene_fam, domlist)]
domlist = domlist[!is.na(domlist)]




#### Plot gene trees ####
annot_from_nn = function(x) {
  
  len = 10
  ida = names(sort(x)[1:len])
  
  # most common annotation
  ann = names(sort(table(stringr::str_split(ida, pattern = "_", simplify = T)[,1]), decreasing = T)[1])
  sum = as.numeric(sort(table(stringr::str_split(ida, pattern = "_", simplify = T)[,1]), decreasing = T)[1] / len)
  
  # closest homolog
  nnn = names(which.min(x))
  nnd = min(x)
  
  # fraction of top hits that are...
  euksu = sum(grepl("^eukoth_", ida)) / len
  arcsu = sum(grepl("^arc_", ida)) / len
  bacsu = sum(grepl("^bac_", ida)) / len
  
  # return
  return(list(ann,sum,nnn,nnd,euksu,arcsu,bacsu))
}


# plot gene trees highlighting taxonomy within archaea
# and orthology of eukaryotic sequences

# record closest group of sequences per sequence? (bacterial, archaeal, or specific eukaryotic OGs)
virt = data.frame()

for (dom in domlist) {
  
  # list phylogenies for this domain
  phyl_li = list.files(pattern = sprintf("euk.%s.*.treefile", dom), path = tre_fo,full.names = T)
  ort_i = ort[grepl(sprintf("^%s\\.", dom),ort$orthogroup),]
  
  pdf(sprintf("results_preeuk/ogs_trees.%s.wtax.pdf", dom), width = 12, height = 12)
  par(mar=c(1.2, 1.2, 1.6, 1.2))
  
  for (phyl_fn in phyl_li) {
    
    print(sprintf("colorful tree %s | %s",dom,phyl_fn))
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
    phyl_nods = data.frame(node = c(phyl$tip.label, phyl$node.label), stringsAsFactors = F)
    phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
    phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)
    
    # taxa in each node
    phyl_nods$taxa = ""
    phyl_nods[phyl_nods$is_tip,"taxa"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "_", simplify = T)[,1]
    # gene name
    phyl_nods$gene = gsub(pattern = "_\\d+-\\d+$", replacement = "",phyl_nods$node)
    
    # add species info
    phyl_nods$gene_id = gsub("^arc_","", phyl_nods$gene)
    phyl_nods$gene_id = gsub("^eukoth_","", phyl_nods$gene_id)
    phyl_nods$gene_id = gsub("^bac_","", phyl_nods$gene_id)
    phyl_nods$gene_id = gsub("_\\d+-\\d+$","", phyl_nods$gene_id)
    # phyl_edge_order = as.character(phyl_edge$ix_edges)
    vtax_i = vtax[vtax$gene %in% phyl_nods$gene_id,]
    phyl_nodt = merge(phyl_nods, vtax_i, by.x = "gene_id", by.y = "gene", all.x = T, all.y=F)
    # add taxonomic info
    taxne_i = taxne[taxne$taxon %in% phyl_nodt$taxon,]
    phyl_nodt = merge(phyl_nodt, taxne_i, by.x = "taxon", by.y = "taxon", all.x = T, all.y = F)
    # reorder
    phyl_nodt = phyl_nodt[match(phyl_nods$node, phyl_nodt$node),]
    
    # get genus data and color string
    phyl_nodt$tax_string = paste(phyl_nodt$superkingdom,phyl_nodt$phylum)
    phyl_nodt$tax_string = sub("Archaea ","", phyl_nodt$tax_string)
    phyl_nodt_genus_tab = table(phyl_nodt$tax_string)
    phyl_nodt_genus_tab = phyl_nodt_genus_tab[names(phyl_nodt_genus_tab) != "NA NA"] 
    phyl_nodt_genus_cla = names(phyl_nodt_genus_tab)
    phyl_nodt_genus_col = viridis::plasma(max(length(phyl_nodt_genus_cla),1), end = 0.90, begin = 0.05)
    
    
    # add orthology
    phyl_nodo = merge(phyl_nods, ort_i, by.x = "gene_id", by.y = "gene", all.x = T, all.y=F, sort = F)
    phyl_nodo = phyl_nodo[match(phyl_nods$node, phyl_nodo$node),]
    phyl_nodo$orthogroup_id = phyl_nodo$orthogroup
    if (sum(!is.na(phyl_nodo$orthogroup)>0)) {
      phyl_nodo_tab = sort(table(phyl_nodo$orthogroup), decreasing=T)
      phyl_nodo_tab_other = sum(phyl_nodo_tab[phyl_nodo_tab < 5])
      phyl_nodo_tab_othen = names(phyl_nodo_tab[phyl_nodo_tab < 5])
      phyl_nodo_tab = phyl_nodo_tab[phyl_nodo_tab >= 5]
      phyl_nodo_tab = c(phyl_nodo_tab, other = phyl_nodo_tab_other)
      phyl_nodo_tab_string = paste(paste(names(phyl_nodo_tab)," n=",phyl_nodo_tab, sep=""), collapse = "\n")
      phyl_nodo_cla = names(phyl_nodo_tab)
      phyl_nodo_col = viridis::viridis(max(length(phyl_nodo_cla),1), end = 1, begin = 0.15)
      phyl_nodo[phyl_nodo$orthogroup %in% phyl_nodo_tab_othen, "orthogroup"] = "other"
    }
    
    # add colors
    phyl_nods$color = NA
    phyl_nods[phyl_nods$taxa %in% c("eukoth"), "color"] = NA
    phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = NA
    phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "slategray4"
    
    # plot unrooted
    plot.phylo(phyl, font=1, type="u", edge.color = "darkgray", root.edge = T, show.tip.label = F, underscore = T)
    tiplabels(col = phyl_nods$color,frame = "none", pch = 1, cex=0.5)
    legend("topright", col=c("blue","deeppink","slategray4"), pch=1, legend = c("eukoth","archaea","bacteria"), cex=0.5)
    title(sprintf("%s\n%s\nn=%i sequences",dom,basename(phyl_fn), length(phyl$tip.label)), cex.main=0.8)
    # plot viral taxonomy nodes
    for (fui in seq_along(phyl_nodt_genus_cla)) {
      gphy_fus_cla_i = phyl_nodt_genus_cla[fui]
      gphy_fus_col_i = phyl_nodt_genus_col[fui]
      tips_w_genus = which(phyl_nodt$tax_string == gphy_fus_cla_i)
      tiplabels(pch = 1, col = gphy_fus_col_i, height = 4, cex = 0.7, tip = tips_w_genus)
      tiplabels(text = paste(gphy_fus_cla_i,phyl_nods$gene[tips_w_genus]),
                tip = tips_w_genus, frame="none", col=alpha(gphy_fus_col_i,0.6), cex=0.4)
    }
    legend("bottomright", col = phyl_nodt_genus_col, legend = paste(names(phyl_nodt_genus_tab),"n =", phyl_nodt_genus_tab), pch=1, cex=0.5)
    # plot eukaryotic orthology nodes
    if (sum(!is.na(phyl_nodo$orthogroup)>0)) {
      for (fui in seq_along(phyl_nodo_cla)) {
        gphy_fus_cla_i = phyl_nodo_cla[fui]
        gphy_fus_col_i = phyl_nodo_col[fui]
        tips_w_og = which(phyl_nodo$orthogroup == gphy_fus_cla_i)
        tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.5, tip = tips_w_og)
      }
      legend("bottomleft", col = phyl_nodo_col, legend = paste(names(phyl_nodo_tab),"n =", phyl_nodo_tab), pch=19, cex=0.5)
    }
    add.scale.bar()
    
    # pairwise distances from phylogeny
    phyl_nodo_tip = phyl_nodo[phyl_nodo$is_tip,]
    gphy_dist = adephylo::distTips(x=phyl, method="patristic")
    gphy_dism = as.matrix(gphy_dist)
    # subset to sequences from OG v. seqs from other homologs
    ix_ogi = which(!grepl("^arc_|^bac_|^eukoth_", rownames(gphy_dism)))
    ix_hom = which( grepl("^arc_|^bac_|^eukoth_", rownames(gphy_dism)))
    # gphy_dist_vir = gphy_dism[ ix_cel, ix_vir ]
    gphy_dist_vir = matrix(
      gphy_dism[ ix_ogi, ix_hom ],
      ncol = length(ix_ogi))
    colnames(gphy_dist_vir) = rownames(gphy_dism)[ix_ogi]
    rownames(gphy_dist_vir) = rownames(gphy_dism)[ix_hom]
    
    # OG-only table
    virt_i = phyl_nodt[!grepl("^arc_|^bac_|^eukoth_", phyl_nodt$gene),c("gene","node")]
    virt_i = virt_i [ virt_i$node %in% phyl$tip.label, ]
    virt_i$core_domain = dom
    virt_i$homology_group = paste(strsplit(basename(phyl_fn), split = "\\.")[[1]][c(2,3,4)], collapse = ".")
    
    # if there are non-viral sequences in the tree...
    if (nrow(gphy_dist_vir)>0) {
      
      # get annotation of closest neighbours in phylogeny
      virt_i_annots = apply(gphy_dist_vir, 2, function(x) annot_from_nn(x) )
      # add to table
      virt_i_annots = data.frame(matrix(unlist(virt_i_annots), nrow=length(virt_i_annots), byrow=TRUE), stringsAsFactors = F)
      colnames(virt_i_annots) = c("closest_group","closest_group_top10frac","closest_node","closest_node_dist","euk_top_frac","arc_top_frac","bac_top_frac")
      virt_i = cbind(virt_i, virt_i_annots)
      
    } else {
      
      # if there are no cellular seqs in the tree, everything is NA
      virt_i$closest_group = NA
      virt_i$closest_group_top10frac = NA
      virt_i$closest_OG_euk = NA
      virt_i$closest_OG_euk_top10frac = NA
      virt_i$closest_node = NA
      virt_i$closest_node_dist = NA
      virt_i$euk_top_frac = NA
      virt_i$arc_top_frac = NA
      virt_i$bac_top_frac = NA
      
    }
    
    # store annotations of each viral sequence
    virt = rbind(virt, virt_i)
    
  }
  
  dev.off()
  
}

# format output table
virt$node = as.character(virt$node)
virt$closest_group_top10frac = as.numeric(as.character(virt$closest_group_top10frac))
virt$euk_top_frac = as.numeric(as.character(virt$euk_top_frac))
virt$arc_top_frac = as.numeric(as.character(virt$arc_top_frac))
virt$bac_top_frac = as.numeric(as.character(virt$bac_top_frac))
virt$closest_group = factor(virt$closest_group, levels = c("eukoth","arc","bac"))

# declare low-freq closest groups as NA
virt$closest_group_filt = factor(virt$closest_group, levels = c("eukoth","unk","arc","bac"))
virt [ virt$bac_top_frac < 0.5 & virt$arc_top_frac < 0.5 & virt$euk_top_frac < 0.5, "closest_group_filt" ] = "unk"

# add virus gene architecture
# virt = merge(virt, arq, all.x = T, all.y = F, by.x = "gene", by.y = "gene", sort = F)






#### Plots per OG ####

# summarise the number of classified genes within each family or OG that have been assigned as 'close' to euks or proks

pdf(file = sprintf("results_preeuk/ogs_closest_homologs.bars.pdf", dom),width = 4, height = 12)
par()$mar
par(mar=c(5.1,8,4.1,2.1))
# globally
# barplots that reflect the number of closest hits of sequences within each gene family
xtab = table(virt$core_domain, virt$closest_group_filt)
xtab_frac = xtab / rowSums(xtab)
xtab_frac = t(xtab_frac)
xtab_col = c("chartreuse3","darkolivegreen2","slategray4","lightgray")
barplot(xtab_frac, horiz = T, xlim=c(0,1), las=1, col=xtab_col, border = "white", cex.axis = 0.6, cex.names = 0.6, ylim = c(0,120))
title(main="all")
legend("topright", cex=0.6, fill = xtab_col, legend = rownames(xtab_frac))

# within each domain
for (dom in domlist) {
  
  # subset
  virt_i = virt [ virt$core_domain == dom, ]
  xtab = table(virt_i$homology_group, virt_i$closest_group_filt)
  xtab_frac = xtab / rowSums(xtab)
  xtab_frac = t(xtab_frac)
  barplot(xtab_frac, horiz = T, xlim=c(0,1), las=1, col=xtab_col, border = "white", cex.axis = 0.6, cex.names = 0.6, ylim = c(0,120))
  title(main=dom)
  legend("topright", cex=0.6, fill = xtab_col, legend = rownames(xtab_frac))
  
}
dev.off()



# same, but heatmaps
# color map
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

pdf(file = sprintf("results_preeuk/ogs_closest_homologs.heatmaps.pdf", dom),width = 3, height = 7)

# within each domain
for (dom in domlist) {
  
  # subset
  virt_i = virt [ virt$core_domain == dom, ]
  xtab = table(virt_i$homology_group, virt_i$closest_group_filt)
  xtab_frac = xtab / rowSums(xtab)

  # heatmap  
  pheatmap(xtab_frac,
           color = col_blue(10), breaks = seq(0,1,length.out = 11),
           cellwidth = 5, cellheight = 5, na_col = "grey",
           cluster_rows = F, cluster_cols = F,
           fontsize = 4,
           main = sprintf("%s", dom),
           border_color = "white", display_numbers = F, number_format = "%.2f")
}
dev.off()


#### Save ####
# save output
write.table(virt, file = "results_preeuk/ogs_summary_preeuk_homologs_annotation.csv", 
            quote = F, sep = "\t", row.names = F)


