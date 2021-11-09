# load libraries
library(viridis)
library(ape)
library(scales)
library(phytools)
library(stringr)
library(data.table)
library(pheatmap)
library(plyr)
# library(MCL)

#### Define input ####

# domains with viral hits
domlist=c("Acetyltransf_1","MOZ_SAS","Hist_deacetyl","SIR2","DOT1","SET","SNF2_N","ASF1_hist_chap","CupinJmjC","GNAT_acetyltr_2","PWWP","MBT","ING","Histone","YEATS","Nucleoplasmin")
# domlist=c("MOZ_SAS","Hist_deacetyl","SIR2","DOT1","SET","ASF1_hist_chap")
# domlist=c("DOT1","SIR2")

# tables
vtax_fn = "../data/seq_Archaea.taxa.csv.gz"    # taxon for each viral seq
taxne_fn = "../data/taxonomy_ranked_wg.tsv.gz" # taxonomy for non-euks
ort_fn = "orthogroups_euk.csv"
cla_fn = "../data/gene_families_hmm.csv"
tre_fo = "results_preeuk/trees_fam/"

# read data
# taxonomy of viral sequences
vtax = fread(cmd=sprintf("zcat %s", vtax_fn), data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
colnames(vtax) = c("gene", "taxon")
# taxonomy of species
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=FALSE, stringsAsFactors=FALSE)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")
taxne$phylum = gsub("Candidatus ", "",taxne$phylum)
# subset to viral sequences for speed
# taxne = taxne[taxne$superkingdom == "Viruses",]
# orthology, gene fam classification, architectures
ort = read.table(ort_fn, sep="\t", header = TRUE, stringsAsFactors = FALSE)
cla = read.table(
  cla_fn, header = FALSE, sep = "\t", 
  col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
  stringsAsFactors = FALSE)


# ordered list of domains
domlist = domlist[match(cla$gene_fam, domlist)]
domlist = domlist[!is.na(domlist)]




#### Plot gene trees ####

annot_from_nn2 = function(x) {
  
  len = 10
  ids = names(sort(x)[1:len])
  ida = ids
  
  # most common annotation
  ann = names(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1])
  sum = as.numeric(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1] / len)
  
  # closest homolog
  nnn = names(which.min(x))
  nnd = min(x)
  
  # fraction of top hits that are...
  arcsu = sum(grepl("^arc_", ida)) / len
  bacsu = sum(grepl("^bac_", ida)) / len
  
  # return
  return(list(ann,sum,nnn,nnd,arcsu,bacsu))
  
}


annot_from_nn3 = function(x) {
  
  len = 10
  ids = names(sort(x)[1:len])
  ida = ids
  ida [ !grepl("^bac|^arc", ida) ] = "eukoth"
  
  
  # most common annotation
  ann = names(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1])
  sum = as.numeric(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1] / len)
  
  # closest homolog
  nnn = names(which.min(x))
  nnd = min(x)
  
  # fraction of top hits that are...
  euksu = sum(ida == "eukoth") / len
  arcsu = sum(grepl("^arc_", ida)) / len
  bacsu = sum(grepl("^bac_", ida)) / len
  
  # return
  return(list(ann,sum,nnn,nnd,euksu,arcsu,bacsu))
  
}


annot_from_nn3_all = function(x) {
  
  len = 10
  ids = names(sort(x)[1:len])
  ida = ids
  ida [ !grepl("^bac|^arc", ida) ] = "eukoth"
  
  
  # most common annotation
  ann = names(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1])
  sum = as.numeric(sort(table(stringr::str_split(ida, pattern = "_", simplify = TRUE)[,1]), decreasing = TRUE)[1] / len)
  
  # closest homolog
  nnn = names(which.min(x))
  nnd = min(x)
  
  # fraction of top hits that are...
  euksu = sum(ida == "eukoth") / len
  arcsu = sum(grepl("^arc_", ida)) / len
  bacsu = sum(grepl("^bac_", ida)) / len
  
  # total number of hits in max support node
  if (ann == "eukoth") {
    annsum = sum(ida == "eukoth")
  } else if (ann == "bac") {
    annsum = sum(grepl("^bac", ida))
  } else if (ann == "arc") {
    annsum = sum(grepl("^arc", ida))
  } else {
    annsum = NA
  }
  
  # return
  return(list(ann,sum,nnn,nnd,euksu,arcsu,bacsu,annsum))
  
}




# plot gene trees highlighting taxonomy within viruses
# and orthology of eukaryotic sequences

# record closest group of sequences per sequence? (bacterial, archaeal, or specific eukaryotic OGs)
virt = data.frame()

for (dom in domlist) {
  
  # list phylogenies for this domain
  phyl_li = list.files(pattern = sprintf("euk.%s.*.treefile", dom), path = tre_fo, full.names = TRUE)
  ort_i = ort[grepl(sprintf("^%s\\.", dom),ort$orthogroup),]
  
  pdf(sprintf("results_preeuk/fam_trees.%s.wtax.pdf", dom), width = 12, height = 12)
  par(mar=c(1.2, 1.2, 1.6, 1.2))
  
  for (phyl_fn in phyl_li) {
    
    # load phylogeny
    phyl = read.tree(phyl_fn)
    phyl = midpoint.root(phyl)
    phyl = ladderize(phyl)
    
    has_euk = any(!grepl( "^bac_|^arc_", phyl$tip.label))
    has_cel = any(grepl( "^bac_|^arc_", phyl$tip.label))
    
    if (has_euk & has_cel) {
      
      print(sprintf("colorful tree %s | %s",dom, phyl_fn))
      
      # dataframe of edges
      phyl_edge           = as.data.frame(phyl$edge)
      colnames(phyl_edge) = c("edge_start","edge_end")
      phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
      phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
      # dataframe of nodes
      phyl_nods = data.frame(node = c(phyl$tip.label, phyl$node.label), stringsAsFactors = FALSE)
      phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
      phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)
      
      # taxa in each node
      phyl_nods$taxa = ""
      phyl_nods[phyl_nods$is_tip,"taxa"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "_", simplify = TRUE)[,1]
      # gene name
      phyl_nods$gene = gsub(pattern = "_\\d+-\\d+$", replacement = "",phyl_nods$node)
      
      # add species info
      phyl_nods$gene_id = gsub("^arc_","", phyl_nods$gene)
      phyl_nods$gene_id = gsub("_\\d+-\\d+$","", phyl_nods$gene_id)
      # phyl_edge_order = as.character(phyl_edge$ix_edges)
      vtax_i = vtax[vtax$gene %in% phyl_nods$gene_id,]
      phyl_nodt = merge(phyl_nods, vtax_i, by.x = "gene_id", by.y = "gene", all.x = TRUE, all.y=FALSE)
      # add taxonomic info
      taxne_i = taxne[taxne$taxon %in% phyl_nodt$taxon,]
      phyl_nodt = merge(phyl_nodt, taxne_i, by.x = "taxon", by.y = "taxon", all.x = TRUE, all.y = FALSE)
      # reorder
      phyl_nodt = phyl_nodt[match(phyl_nods$node, phyl_nodt$node),]
      
      # get genus data and color string
      phyl_nodt$tax_string = paste(phyl_nodt$superkingdom,phyl_nodt$phylum)
      phyl_nodt$tax_string = sub("Archaea ","", phyl_nodt$tax_string)
      phyl_nodt_genus_tab = table(phyl_nodt$tax_string)
      phyl_nodt_genus_tab = phyl_nodt_genus_tab[names(phyl_nodt_genus_tab) != "NA NA"] 
      phyl_nodt_genus_cla = names(phyl_nodt_genus_tab)
      phyl_nodt_genus_col = viridis::plasma(max(length(phyl_nodt_genus_cla),1), end = 0.90, begin = 0.2)
      
      
      # add orthology
      phyl_nodo = merge(phyl_nods, ort_i, by.x = "gene", by.y = "gene", all.x = TRUE, all.y=FALSE, sort = FALSE)
      phyl_nodo = phyl_nodo[match(phyl_nods$node, phyl_nodo$node),]
      phyl_nodo$orthogroup_id = phyl_nodo$orthogroup
      if (sum(!is.na(phyl_nodo$orthogroup) > 0)) {
        phyl_nodo_tab = sort(table(phyl_nodo$orthogroup), decreasing=TRUE)
        phyl_nodo_tab_other = sum(phyl_nodo_tab[phyl_nodo_tab < 10])
        phyl_nodo_tab_othen = names(phyl_nodo_tab[phyl_nodo_tab < 10])
        phyl_nodo_tab = phyl_nodo_tab[phyl_nodo_tab >= 10]
        phyl_nodo_tab = c(phyl_nodo_tab, other = phyl_nodo_tab_other)
        phyl_nodo_tab_string = paste(paste(names(phyl_nodo_tab)," n=",phyl_nodo_tab, sep=""), collapse = "\n")
        phyl_nodo_cla = names(phyl_nodo_tab)
        phyl_nodo_col = viridis::viridis(max(length(phyl_nodo_cla),1), end = 1, begin = 0.15)
        phyl_nodo[phyl_nodo$orthogroup %in% phyl_nodo_tab_othen, "orthogroup"] = "other"
      }
      
      
      # add colors
      phyl_nods$color = NA
      phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = NA
      phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "slategray4"
      phyl_nods[phyl_nods$taxa %in% c("vir"), "color"] = NA
      
      
      # plot unrooted
      plot.phylo(phyl, font=1, type="u", edge.color = "darkgray", root.edge = TRUE, show.tip.label = FALSE, underscore = TRUE)
      tiplabels(col = phyl_nods$color,frame = "none", pch = 1, cex=0.5)
      legend("topright", col=c("deeppink","blue","slategray4"), pch=1, legend = c("eukaryotes","archaea","bacteria"), cex=0.5)
      title(sprintf("%s\n%s\nn=%i sequences",dom,basename(phyl_fn), length(phyl$tip.label)), cex.main=0.8)
      # plot archaeal taxonomy nodes
      for (fui in seq_along(phyl_nodt_genus_cla)) {
        gphy_fus_cla_i = phyl_nodt_genus_cla[fui]
        gphy_fus_col_i = phyl_nodt_genus_col[fui]
        tips_w_genus = which(phyl_nodt$tax_string == gphy_fus_cla_i)
        tiplabels(pch = 1, col = gphy_fus_col_i, height = 4, cex = 0.5, tip = tips_w_genus)
        tiplabels(text = paste(gphy_fus_cla_i,phyl_nods$gene[tips_w_genus]),
                  tip = tips_w_genus, frame="none", col=alpha(gphy_fus_col_i,0.6), cex=0.4)
      }
      legend("bottomright", col = phyl_nodt_genus_col, legend = paste(names(phyl_nodt_genus_tab),"n =", phyl_nodt_genus_tab), pch=1, cex=0.5)
      # plot eukaryotic orthology nodes
      if (sum(!is.na(phyl_nodo$orthogroup) > 0)) {
        for (fui in seq_along(phyl_nodo_cla)) {
          gphy_fus_cla_i = phyl_nodo_cla[fui]
          gphy_fus_col_i = phyl_nodo_col[fui]
          tips_w_og = which(phyl_nodo$orthogroup == gphy_fus_cla_i)
          tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.5, tip = tips_w_og)
          if (gphy_fus_cla_i != "other") {
            tiplabels(text = phyl_nods$gene[tips_w_og],
                      tip = tips_w_og, frame="none", col=alpha(gphy_fus_col_i,0.6), cex=0.4)
          }
        }
        legend("bottomleft", col = phyl_nodo_col, legend = paste(names(phyl_nodo_tab),"n =", phyl_nodo_tab), pch=19, cex=0.5)
      }
      add.scale.bar()
      
      # add colors
      phyl_nods$color = "deeppink"
      phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = "blue"
      phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "slategray4"
      
      ape::plot.phylo(phyl, font=1, type="u", edge.color = "darkgray", root.edge = TRUE, show.tip.label = TRUE, underscore = TRUE,
                 show.node.label = TRUE, cex=0.2, tip.color = alpha(phyl_nods$color,0.5))
      tiplabels(pch=19, cex=0.1, col=phyl_nods$color)
      add.scale.bar()
      
      
      
      ### PREP: PAIRWISE DISTANCES
      # pairwise distances from phylogeny
      phyl_nodo_tip = phyl_nodo[phyl_nodo$is_tip,]
      # pairwise patritic distance:
      gphy_dism = ape::cophenetic.phylo(phyl)
      
      # euk-only table
      virt_i = phyl_nodt[
        !grepl("^bac_|^arc_", phyl_nodt$gene) & phyl_nodt$node %in% phyl$tip.label,
        c("gene","node")]
      virt_i$core_domain = dom
      virt_i$homology_group = paste(strsplit(basename(phyl_fn), split = "\\.")[[1]][c(2,4)], collapse = ".")
      virt_i = merge(virt_i, ort_i, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
      
      
      # annotate each gene to its closest internal cluster OR to prokaryotes
      list_internal_clu = unique(virt_i$orthogroup)
      list_internal_clu = list_internal_clu[!is.na(list_internal_clu)]
      list_internal_clu_df = data.frame()
      if (length(list_internal_clu) > 0) {
        for (cli in list_internal_clu) {
          
          ix_ref = which(
            !grepl("^bac_|^arc_", rownames(gphy_dism))
            & rownames(gphy_dism) %in% virt_i[virt_i$orthogroup == cli,]$node 
          )
          ix_que = which(
            grepl("^bac_|^arc_", rownames(gphy_dism))
            | !rownames(gphy_dism) %in% virt_i[virt_i$orthogroup == cli,]$node 
          )
          gphy_dist_que = matrix(
            gphy_dism[ ix_que, ix_ref ],
            ncol = length(ix_ref))
          colnames(gphy_dist_que) = rownames(gphy_dism)[ix_ref]
          rownames(gphy_dist_que) = rownames(gphy_dism)[ix_que]
          
          virt_i_annots = apply(gphy_dist_que, 2, function(x) annot_from_nn3_all(x) )
          virt_i_annotd = data.frame(matrix(unlist(virt_i_annots), nrow=length(virt_i_annots), byrow=TRUE), stringsAsFactors = FALSE)
          colnames(virt_i_annotd) = c("closest_group","closest_group_top10frac","closest_node","closest_node_dist","euk_top_frac","arc_top_frac","bac_top_frac","closest_node_num")
          virt_i_annotd$node = names(virt_i_annots)
          
          # store dataframe
          list_internal_clu_df = rbind(list_internal_clu_df, virt_i_annotd)
        }
        virt_i = merge(virt_i, list_internal_clu_df, by.x = "node", by.y = "node", all.x = TRUE, all.y = FALSE, sort = FALSE)
      }
      
      ### FINAL
      # store annotations
      virt = plyr::rbind.fill(virt, virt_i)
      
      
    } # closes IF that checks if this individual tree has euk and prok genes
  } # closes FOR LOOP of trees within this family
  
  dev.off()
  
}

# format output table
virt$node = as.character(virt$node)
virt$arc_top_frac = as.numeric(as.character(virt$arc_top_frac))
virt$bac_top_frac = as.numeric(as.character(virt$bac_top_frac))
virt$euk_top_frac = as.numeric(as.character(virt$euk_top_frac))
virt$closest_node_num =  as.numeric(as.character(virt$closest_node_num))
virt$arc_top_frac [ is.na(virt$arc_top_frac) ] = 0
virt$bac_top_frac [ is.na(virt$bac_top_frac) ] = 0
virt$euk_top_frac [ is.na(virt$euk_top_frac) ] = 0
virt$closest_group_factor = factor(virt$closest_group, levels = c("eukoth","arc","bac","unk"))
virt [ is.na(virt$closest_group_factor) , "closest_group_factor" ] = "unk"
virt [ is.na(virt$closest_node_num) , "closest_node_num" ] = 0

# declare low-freq closest groups as NA
virt$closest_group_factor_filt = factor(virt$closest_group_factor, levels = c("eukoth","unk","arc","bac"))
virt [ virt$bac_top_frac < 0.5 & virt$arc_top_frac < 0.5 & virt$euk_top_frac < 0.5, "closest_group_factor_filt" ] = "unk"
virt [ virt$closest_node_num <= 5, "closest_group_factor_filt" ] = "unk"


# core domain factors
virt$core_domain = factor(virt$core_domain, levels = cla$gene_fam)
virt$core_domain = droplevels(virt$core_domain)


#### Plots per OG ####

# summarise the number of classified genes within each family or OG that have been assigned as 'close' to euks or proks
pdf(file = sprintf("results_preeuk/summary_closest_homologs.bars.pdf", dom),width = 6, height = 12)
par(mar=c(5.1,22,4.1,2.1))
# globally
# barplots that reflect the number of closest hits of sequences within each gene family
xtab = table(virt$core_domain, virt$closest_group_factor_filt)
xtab_frac = xtab / rowSums(xtab)
xtab_frac = t(xtab_frac)
xtab_col = c("chartreuse3","darkolivegreen2","slategray4","lightgray")
barplot(xtab_frac, horiz = TRUE, xlim=c(0,1), las=1, col=xtab_col, border = "white", cex.axis = 0.6, cex.names = 0.6, ylim = c(0,100))
title(main="all")
legend("topright", cex=0.6, fill = xtab_col, legend = rownames(xtab_frac))

# ignore rare OGs
ort_tab = sort(table(ort$orthogroup), decreasing = TRUE)
ort_tab = ort_tab [ ort_tab >= 10 ]

# within each gene family
for (dom in domlist) {
  
  # subset
  virt_i = virt [ virt$core_domain == dom, ]
  ort_in_i = names(ort_tab)[ names(ort_tab) %in% virt_i$orthogroup ]
  virt_i = virt_i [ virt_i$orthogroup %in% ort_in_i, ]
  
  # cross-tabulate
  xtab = table(virt_i$orthogroup, virt_i$closest_group_factor_filt)
  xtab_frac = xtab
  xtab_frac = t(xtab_frac)
  xtab_frac = xtab_frac[ ,order(colSums(xtab_frac), decreasing = TRUE)  ]
  barplot(xtab_frac, horiz = TRUE, xlim=c(0,500), las=1, col=xtab_col, border = "white", cex.axis = 0.6, cex.names = 0.6, ylim = c(0,100))
  title(main=dom)
  legend("topright", cex=0.6, fill = xtab_col, legend = rownames(xtab_frac))
  
}
dev.off()





#### Save ####

# save output
write.table(virt, file = "results_preeuk/summary_preeuk_families.annotation_per_gene.csv", 
            quote = FALSE, sep = "\t", row.names = FALSE)


# save per-OG table
virt_sum = virt [ , c("orthogroup", "closest_group", "closest_group_top10frac") ]
virt_sum$closest_group_top10frac = as.numeric(virt_sum$closest_group_top10frac)
virt_agg = aggregate(closest_group_top10frac ~ orthogroup + closest_group, data=virt_sum, sum)
virt_agg = virt_agg [ order(virt_agg$orthogroup, virt_agg$closest_group_top10frac, decreasing = TRUE) , ]

# aggregate: closest outgroup
virt_agg_top = virt_agg[ !duplicated(virt_agg$orthogroup) , ]
colnames(virt_agg_top)[3] = "closest_group_total_frac"

# aggregate: second-closest outgroup
virt_agg_oth = virt_agg[  duplicated(virt_agg$orthogroup) , ]
virt_agg_oth = virt_agg_oth[ !duplicated(virt_agg_oth$orthogroup) , ]
colnames(virt_agg_oth)[2:3] = c("second_group","second_group_total_frac")

# merge
virt_agg_all = merge(virt_agg_top, virt_agg_oth, by = "orthogroup", all.x = TRUE, all.y = FALSE)
virt_agg_all$closest_score = virt_agg_all$closest_group_total_frac / (virt_agg_all$second_group_total_frac + virt_agg_all$closest_group_total_frac)
virt_agg_all [ is.na(virt_agg_all$closest_score) , "closest_score" ] = 1
virt_agg_all = virt_agg_all [ , c("orthogroup", "closest_group", "closest_score")]
virt_agg_all$orthogroup_key = gsub(":.*","", virt_agg_all$orthogroup)

# add info about presence probability at LECA
euk_prob = read.table("orthogroups_euk.ancestral.posteriors_pres.csv")
euk_prob_leca = data.frame(orthogroup_key = rownames(euk_prob), LECA_probability = euk_prob$Eukaryota)
is_null_probability = is.na(euk_prob_leca$LECA_probability)
# if presence probability in LECA not available, retrieve from Wagner parsimony
wag_prob = read.table("orthogroups_euk.ancestral.wagner_g5_pres.csv", comment.char = "", sep = "\t", header = TRUE)
wag_prob[,1] = NULL
wag_prob$name = gsub(":.*","", wag_prob$name)
rownames(wag_prob) = wag_prob$name
wag_prob = wag_prob [ euk_prob_leca$orthogroup_key, ]
euk_prob_leca$LECA_probability [ is_null_probability ] = wag_prob$Eukaryota [ is_null_probability ] > 0
euk_prob_leca$is_probability = !is_null_probability
virt_agg_tot = merge(virt_agg_all, euk_prob_leca, by = "orthogroup_key", all.x = TRUE, all.y = TRUE)


# ort dict
ord = data.frame(orthogroup = ort$orthogroup, orthogroup_key = gsub(":.*","",ort$orthogroup))
ord = ord[!duplicated(ord) , ]
virt_agg_tot = virt_agg_tot [ , c("orthogroup_key", "closest_group", "closest_score", "LECA_probability", "is_probability")]
virt_agg_tot = merge(virt_agg_tot, ord, by = "orthogroup_key", all.x = TRUE, all.y = FALSE)
virt_agg_tot = virt_agg_tot [ , c("orthogroup", "closest_group", "closest_score", "LECA_probability", "is_probability")]

# save
write.table(virt_agg_tot, file = "results_preeuk/summary_preeuk_families.annotation_per_OG.csv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

message("all done")