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
domlist=c("SET","DOT1","SIR2","Hist_deacetyl","MOZ_SAS","GNAT_acetyltr_2","CupinJmjC","SNF2_N")

# tables
vtax_fn = "../data/seq_Archaea.taxa.csv.gz"    # taxon for each viral seq
taxne_fn = "../data/taxonomy_ranked_wg.tsv.gz" # taxonomy for non-euks
ort_fn = "orthogroups_euk.csv"
cla_fn = "../data/gene_families_hmm.csv"
arq_fn = "gene_sequences/noneuk_arc.pfamscan_archs.csv"

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

# first, plot trees highlighting viral seqs
for (dom in domlist) {
  
  print(sprintf("simple tree %s",dom))
  
  # list phylogenies for this domain
  phyl_li = list.files(pattern = sprintf("arc.%s.*.treefile", dom), path = "results_archaea/trees/",full.names = T)
  ort_i = ort[grepl(sprintf("^%s\\.", dom),ort$orthogroup),]
  
  pdf(sprintf("results_archaea/trees.%s.pdf", dom), width = 12, height = 12)
  par(mar=c(1.2, 1.2, 1.6, 1.2))
  
  for (phyl_fn in phyl_li) {
    
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
    
    # add orthology
    phyl_nodo = merge(phyl_nods, ort_i, by.x = "gene", by.y = "gene", all.x = F)
    if (nrow(phyl_nodo)>0) {
      phyl_nodo_tab = sort(table(phyl_nodo$orthogroup), decreasing=T)
      phyl_nodo_tab = phyl_nodo_tab[phyl_nodo_tab>10]
      phyl_nodo_tab_string = paste(paste(names(phyl_nodo_tab)," n=",phyl_nodo_tab, sep=""), collapse = "\n")
    }
    
    # add colors
    phyl_nods$color = "deeppink"
    phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = "blue"
    phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "slategray4"
    
    # plot unrooted
    plot.phylo(phyl, font=1, type="u", edge.color = "darkgray", root.edge = T, show.tip.label = F, underscore = T)
    tiplabels(col = phyl_nods$color,frame = "none", pch = 19, cex=0.7)
    tiplabels(text = phyl_nods[phyl_nods$taxa=="arc","gene"],tip = which(phyl_nods$taxa == "arc"), frame="none", col=alpha("blue",0.6), cex=0.4)
    add.scale.bar()
    title(main = sprintf("%s\n%s\nn=%i sequences",dom,basename(phyl_fn), length(phyl$tip.label)))
    if (nrow(phyl_nodo)>0) {
      title(sub = phyl_nodo_tab_string, cex.sub=0.5)
    }
    legend("topright", col=c("deeppink","blue","slategray4"), pch=19, legend = c("eukaryotes","archaea","bacteria"), cex=0.5)
    
  }
  
  dev.off()
  
}

annot_from_nn = function(x) {
  
  len = min(10,length(x))
  ids = names(sort(x)[1:len])
  ida = ids
  ida [ ida %in% as.character(phyl_nodo_tip$node[!is.na(phyl_nodo_tip$orthogroup)]) ] = "euk"
  
  # most common annotation
  ann = names(sort(table(stringr::str_split(ida, pattern = "_", simplify = T)[,1]), decreasing = T)[1])
  sum = as.numeric(sort(table(stringr::str_split(ida, pattern = "_", simplify = T)[,1]), decreasing = T)[1] / len)
  
  # most common orthogroup
  ogt = phyl_nodo[phyl_nodo$node %in% ids , "orthogroup_id"]
  oga = names(sort(table(ogt), decreasing = T)[1])
  ogs = as.numeric(sort(table(ogt), decreasing = T)[1] / len)
  if (is.null(oga)) { 
    oga = NA
    ogs = NA
  }
  
  # closest homolog
  nnn = names(which.min(x))
  nnd = min(x)
  
  # fraction of top hits that are...
  euksu = sum(ida == "euk") / len
  bacsu = sum(grepl("^bac_", ida)) / len
  
  # return
  return(list(ann,sum,oga,ogs,nnn,nnd,euksu,bacsu))
}



# second, plot gene trees highlighting taxonomy within viruses
# and orthology of eukaryotic sequences

# record closest group of sequences per sequence? (bacterial, archaeal, or specific eukaryotic OGs)
virt = data.frame()

for (dom in domlist) {
  
  print(sprintf("colorful tree %s",dom))
  
  # list phylogenies for this domain
  phyl_li = list.files(pattern = sprintf("arc.%s.*.treefile", dom), path = "results_archaea/trees/",full.names = T)
  ort_i = ort[grepl(sprintf("^%s\\.", dom),ort$orthogroup),]
  
  pdf(sprintf("results_archaea/trees.%s.wtax.pdf", dom), width = 12, height = 12)
  par(mar=c(1.2, 1.2, 1.6, 1.2))
  
  for (phyl_fn in phyl_li) {
    
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
    phyl_nodo = merge(phyl_nods, ort_i, by.x = "gene", by.y = "gene", all.x = T, all.y=F, sort = F)
    phyl_nodo = phyl_nodo[match(phyl_nods$node, phyl_nodo$node),]
    phyl_nodo$orthogroup_id = phyl_nodo$orthogroup
    if (sum(!is.na(phyl_nodo$orthogroup)>0)) {
      phyl_nodo_tab = sort(table(phyl_nodo$orthogroup), decreasing=T)
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
    phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = "blue"
    phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "slategray4"
    phyl_nods[phyl_nods$taxa %in% c("vir"), "color"] = NA
    
    # plot unrooted
    plot.phylo(phyl, font=1, type="u", edge.color = "darkgray", root.edge = T, show.tip.label = F, underscore = T)
    tiplabels(col = phyl_nods$color,frame = "none", pch = 1, cex=0.5)
    legend("topright", col=c("deeppink","blue","slategray4"), pch=1, legend = c("eukaryotes","archaea","bacteria"), cex=0.5)
    title(sprintf("%s\n%s\nn=%i sequences",dom,basename(phyl_fn), length(phyl$tip.label)), cex.main=0.8)
    # plot viral taxonomy nodes
    for (fui in seq_along(phyl_nodt_genus_cla)) {
      gphy_fus_cla_i = phyl_nodt_genus_cla[fui]
      gphy_fus_col_i = phyl_nodt_genus_col[fui]
      tips_w_genus = which(phyl_nodt$tax_string == gphy_fus_cla_i)
      tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.7, tip = tips_w_genus)
      tiplabels(text = paste(gphy_fus_cla_i,phyl_nods$gene[tips_w_genus]),
                tip = tips_w_genus, frame="none", col=alpha(gphy_fus_col_i,0.6), cex=0.4)
    }
    legend("bottomright", col = phyl_nodt_genus_col, legend = paste(names(phyl_nodt_genus_tab),"n =", phyl_nodt_genus_tab), pch=19, cex=0.5)
    # plot eukaryotic orthology nodes
    if (sum(!is.na(phyl_nodo$orthogroup)>0)) {
      for (fui in seq_along(phyl_nodo_cla)) {
        gphy_fus_cla_i = phyl_nodo_cla[fui]
        gphy_fus_col_i = phyl_nodo_col[fui]
        tips_w_og = which(phyl_nodo$orthogroup == gphy_fus_cla_i)
        tiplabels(pch = 1, col = gphy_fus_col_i, height = 4, cex = 0.5, tip = tips_w_og)
      }
      legend("bottomleft", col = phyl_nodo_col, legend = paste(names(phyl_nodo_tab),"n =", phyl_nodo_tab), pch=1, cex=0.5)
    }
    add.scale.bar()
    
    # pairwise distances from phylogeny
    phyl_nodo_tip = phyl_nodo[phyl_nodo$is_tip,]
    # phyl$edge.length = rep(1,length(phyl$edge.length))
    # gphy_dist = ape::cophenetic.phylo(phyl)
    gphy_dist = adephylo::distTips(x=phyl, method="patristic")
    gphy_dism = as.matrix(gphy_dist)
    # subset to virus v. cellular sequences
    ix_vir = which(grepl("^arc_", rownames(gphy_dism)))
    ix_cel = which(!grepl("^arc_", rownames(gphy_dism)))
    # gphy_dist_vir = gphy_dism[ ix_cel, ix_vir ]
    gphy_dist_vir = matrix(
      gphy_dism[ ix_cel, ix_vir ],
      ncol = length(ix_vir))
    colnames(gphy_dist_vir) = rownames(gphy_dism)[ix_vir]
    rownames(gphy_dist_vir) = rownames(gphy_dism)[ix_cel]
    
    # gphy_dist_vir = raw_pairwise_dist(
    #   phy = phyl, 
    #   tip.x = phyl$tip.label[!grepl("^vir_",phyl$tip.label)], 
    #   tip.y = phyl$tip.label[grepl("^vir_",phyl$tip.label)])
    
    # virus-only table
    virt_i = phyl_nodt[grepl("^arc_", phyl_nodt$gene),c("gene","node","taxon","genus","phylum","kingdom")]
    virt_i$core_domain = dom
    virt_i$homology_group = paste(strsplit(basename(phyl_fn), split = "\\.")[[1]][c(2,4)], collapse = ".")
    
    # if there are non-viral sequences in the tree...
    if (nrow(gphy_dist_vir)>0) {
      
      # get annotation of closest neighbours in phylogeny
      virt_i_annots = apply(gphy_dist_vir, 2, function(x) annot_from_nn(x) )
      # add to table
      virt_i_annots = data.frame(matrix(unlist(virt_i_annots), nrow=length(virt_i_annots), byrow=TRUE), stringsAsFactors = F)
      colnames(virt_i_annots) = c("closest_group","closest_group_top10frac","closest_OG_euk","closest_OG_euk_top10frac","closest_node","closest_node_dist","euk_top_frac","bac_top_frac")
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
virt$closest_OG_euk_top10frac = as.numeric(as.character(virt$closest_OG_euk_top10frac))
virt$euk_top_frac = as.numeric(as.character(virt$euk_top_frac))
virt$bac_top_frac = as.numeric(as.character(virt$bac_top_frac))

# add virus gene architecture
virt = merge(virt, arq, all.x = T, all.y = F, by.x = "gene", by.y = "gene", sort = F)

# save output
write.table(virt, file = "results_archaea/summary_archaeal_homologs_annotation.csv", 
            quote = F, sep = "\t", row.names = F)


