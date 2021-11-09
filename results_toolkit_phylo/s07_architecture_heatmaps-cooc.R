# load libraries
library(scales)
library(stringr)
library(pheatmap)
library(ape)
library(tidyr)

graphics.off()

#### Define input ####

clas_fn = "../data/gene_families_hmm.csv"
arqs_fo = "results_domains/"
phyl_fn = "data/species_tree.newick"
tedo_fn = "../data/transposable_element_domains.csv"
tree_fo = "gene_trees//"
garq_fn = "architectures.csv"
gcla_fn = "gene_counts/euk_genecounts.csv"
taxe_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"   # taxonomic info for eukaryotes


col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

# load tables
cla = read.table(clas_fn, header = F, sep = "\t", 
                 col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
                 stringsAsFactors = F)

# species tree for species ordering
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label

# load taxonomy and find out where to add spacing in heatmaps
taxe = read.table(taxe_fn, sep = "\t", header = T, stringsAsFactors = T)
taxe = taxe[match(sps_order, taxe$Species),]
sps_gap_ixs = c(1,1+which(diff(as.numeric(taxe$Group))!=0)) - 1


# list of TEs
te_doms = read.table(tedo_fn, header = T, stringsAsFactors = F, sep = "\t")
list_te_doms = te_doms$domain

# regexp of TE domains
list_te_doms_rexp = paste(list_te_doms, collapse = "\\b|\\b")

#### Plot heatmaps per family ####

list_fams = cla$gene_fam
list_fams_dom_regexp = lapply(stringr::str_split(cla$domains, ",", simplify = F), function(x) paste(x, collapse = "\\b|\\b"))
names(list_fams_dom_regexp) = list_fams

for (fam in list_fams) {
  
  # load arqs for that gene family
  if (file.exists(sprintf("%s/net.%s.network_genes.csv",arqs_fo, fam))) {
    arq = read.table(sprintf("%s/net.%s.network_genes.csv",arqs_fo, fam), sep = "\t", stringsAsFactors = F, header = T)
    arq = arq[arq$family == fam,]
    
    if (nrow(arq)>0) {
      
      print(sprintf("architectures %s", fam))
      
      # add species
      arq$species = stringr::str_split(arq$gene, pattern = "_", simplify = T)[,1]
      arq$species = factor(arq$species, levels = sps_order)
      
      # simplify architectures (collapse **consecutive** repeats)
      simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
      simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
      arq$simarq = simple_arqs_vect
      
      #### Architectures per family ####
      # count architectures per species, for the whole gene family
      arm = data.frame(xtabs(formula = ~ simarq + species, data=arq))
      arm = tidyr::spread(arm, species, Freq)
      rownames(arm) = arm[,1]
      arm = arm [,-1]
      arm_gene_counts = colSums(arm)
      arm = rbind(arm_gene_counts, arm)
      rownames(arm)[1] = "Gene presence"
      arm_rowlabs = stringr::str_trunc(rownames(arm), width = 60, ellipsis = "...")
      arm_rowlabs = stringr::str_pad(arm_rowlabs, width = 60, side = "right")
      # order architectures
      arm_ord = order(rowSums(arm>0) , apply(arm, 1, which.max), decreasing = T)
      
      # which domains are TE domains?
      arm_wte_bool = grepl(list_te_doms_rexp, rownames(arm))
      
      # # plot architectures heatmap
      # pdf(file = sprintf("%s/table.%s.dom_arqs_per_fam.pdf", arqs_fo, fam),width = 16, height = max(4,nrow(arm)/10))
      # pheatmap(arm[arm_ord,],  
      #          color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
      #          gaps_row = 1,
      #          cellwidth = 4, cellheight = 4, na_col = "grey", 
      #          cluster_rows = F, cluster_cols = F,
      #          fontsize = 5, 
      #          annotation_row = data.frame(is_TE = arm_wte_bool*1, row.names = rownames(arm)),
      #          gaps_col = sps_gap_ixs, 
      #          labels_row = arm_rowlabs[arm_ord],
      #          main = sprintf("%s architectures", fam),
      #          border_color = "white", display_numbers = F, number_format = "%i")
      # 
      # # focus on architectures with TE domains
      # if (nrow(arm[arm_wte_bool,]) > 0) {
      #   arm_wte = arm[arm_wte_bool,]
      #   arm_wte_ord = order(rowSums(arm_wte>0),decreasing = T)
      #   pheatmap(arm_wte[arm_wte_ord,],  
      #            color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
      #            cellwidth = 4, cellheight = 4, na_col = "grey", 
      #            cluster_rows = F, cluster_cols = F,
      #            fontsize = 5, 
      #            labels_row = arm_rowlabs[arm_wte_bool][arm_wte_ord],
      #            main = sprintf("%s architectures with TE domains", fam),
      #            border_color = "white", display_numbers = F, number_format = "%i")
      # }
      # dev.off()
      
      #### Domains per family ####
      # obtain domain duets
      # simple_arqs_list_duets = sapply(simple_arqs_list,
      #                                 FUN = function(x) {
      #                                   u=unique(sort(x))
      #                                   ci=gtools::combinations(n=length(u), r=2, repeats.allowed = T)
      #                                   co=sort(apply(ci, 1, FUN = function(c) paste(u[c], collapse = ",")))
      #                                 })
      # or rather, use simply domain presence?
      simple_arqs_list_duets = sapply(simple_arqs_list, FUN = function(x) unique(x) )
      names(simple_arqs_list_duets) = arq$gene
      ard = stack(simple_arqs_list_duets)
      colnames(ard) = c("domains","gene")
      # add species & orthogroup
      ard = merge(ard, arq, by="gene", all.y = F)
      
      # regexp of family-defining domains
      list_fam_doms_rexp = paste(str_split(cla[cla$gene_fam == fam,"domains"], ",", simplify = T)[1,], collapse = "|")
      # subset to duets involving TE domains and family-defining domains
      ard_wte = ard[grepl(list_fams_dom_regexp[fam], ard$domains),]
      ard_wte = ard_wte[grepl(list_fams_dom_regexp[fam], ard_wte$domains),]
      
      # arm = data.frame(xtabs(formula = ~ domains + species, data=ard_wte))
      # if (nrow(arm)>0) {
      #   arm = tidyr::spread(arm, species, Freq)
      #   rownames(arm) = arm[,1]
      #   arm = arm [,-1]
      #   arm_rowlabs = stringr::str_trunc(rownames(arm), width = 60, ellipsis = "...")
      #   arm_rowlabs = stringr::str_pad(arm_rowlabs, width = 60, side = "right")
      #   # order architectures
      #   if (length(arm) > 2) {
      #     arm_ord = order(rowSums(arm>0) , apply(arm, 1, which.max), decreasing = T)
      #     arp = arm[arm_ord,]
      #   } else {
      #     arm_ord = 1:length(arm)
      #     arp = arm[arm_ord]
      #   }
      #   
      #   # plot architectures heatmap: duets only
      #   # first, which duets have a TE domain?
      #   arp_wte_bool = grepl(list_te_doms_rexp, rownames(arp))
      #   pdf(file = sprintf("%s/table.%s.dom_per_fam.pdf", arqs_fo, fam),width = 16, height = max(4,nrow(arp)/10))
      #   pheatmap(arp,  
      #            color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
      #            cellwidth = 4, cellheight = 4, na_col = "grey", 
      #            cluster_rows = F, cluster_cols = F,
      #            fontsize = 5, number_color = "red", fontsize_number = 3, 
      #            gaps_col = sps_gap_ixs, 
      #            annotation_row = data.frame(is_TE = arp_wte_bool*1, row.names = rownames(arp)),
      #            labels_row = arm_rowlabs[arm_ord],
      #            labels_col = stringr::str_trunc(colnames(arm), width = 60, ellipsis = "..."),
      #            main = sprintf("%s domain duets", fam),
      #            border_color = "white", display_numbers = F, number_format = "%i")
      #   dev.off()
      # }
      
      #### Duets per OG ####
      list_ogs = names(sort(table(arq$orthogroup), decreasing = T))
      arm_cat = data.frame()
      arm_cat_breaks = c()
      pdf(file = sprintf("%s/table.%s.dom_per_OG.pdf", arqs_fo, fam),width = 14, height = 6)
      for (ogi in list_ogs) {
        
        # count architectures per species
        ari = ard[ ard$orthogroup == ogi , ]
        arm = data.frame(xtabs(formula = ~ domains + species, data=ari))
        arm = tidyr::spread(arm, species, Freq)
        rownames(arm) = arm[,1]
        arm = arm [,-1]
        # order architectures
        arm_gene_counts = colSums(arm)
        arm = rbind(arm_gene_counts, arm)
        rownames(arm)[1] = "Gene presence"
        arm_rowlabs = stringr::str_trunc(rownames(arm), width = 60, ellipsis = "...")
        arm_rowlabs = stringr::str_pad(arm_rowlabs, width = 60, side = "right")
        arm_ord = order(rowSums(arm>0) , apply(arm, 1, which.max), decreasing = T)
        
        # plot duets heatmap
        arp_wte_bool = grepl(list_te_doms_rexp, rownames(arm[arm_ord,]))
        pheatmap(arm[arm_ord,], 
                 color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
                 cellwidth = 4, cellheight = 4, na_col = "grey", gaps_row = 1,
                 cluster_rows = F, cluster_cols = F,
                 fontsize = 5,
                 annotation_row = data.frame(is_TE = arp_wte_bool*1, row.names = rownames(arm[arm_ord,])),
                 gaps_col = sps_gap_ixs, 
                 labels_row = arm_rowlabs[arm_ord],
                 main = sprintf("%s duets, OG %s", fam, ogi),
                 border_color = "white", display_numbers = F, number_format = "%i")
        
        # if there is a TE in this OG, store it for a later concatenated table
        if (any(arp_wte_bool)) {
          
          arm_f = arm[arm_ord,]
          # arm_f = arm_f[rownames(arm_f) != "Gene presence",]
          rownames(arm_f) = paste(
            stringr::str_trunc(ogi, width = 40, ellipsis = "..."),  
            rownames(arm_f), sep = " | ")
          arm_cat = rbind(arm_cat, arm_f)
          arm_cat_breaks = c(arm_cat_breaks, nrow(arm_cat))
          
        }
      }
      dev.off()    
      
      #### Duets per OG: concatenate TE tables ####
      if (nrow(arm_cat) > 0) {
        pdf(file = sprintf("%s/table.%s.dom_per_OG-with_TE.pdf", arqs_fo, fam),width = 16, height = max(4,nrow(arm_cat)/10))
        arp_wte_bool = grepl(list_te_doms_rexp, rownames(arm_cat))
        pheatmap(arm_cat, 
                 color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
                 cellwidth = 4, cellheight = 4, na_col = "grey", 
                 cluster_rows = F, cluster_cols = F,
                 fontsize = 5,
                 annotation_row = data.frame(is_TE = arp_wte_bool*1, row.names = rownames(arm_cat)),
                 gaps_col = sps_gap_ixs, gaps_row = arm_cat_breaks,
                 main = sprintf("%s duets, OGs with TE insertions", fam),
                 border_color = "white", display_numbers = F, number_format = "%i")
        dev.off()
      }
      
      #### Duets per OG ####
      list_ogs = names(sort(table(arq$orthogroup), decreasing = T))
      arm_cat = data.frame()
      arm_cat_breaks = c()
      pdf(file = sprintf("%s/table.%s.dom_arqs_per_OG.pdf", arqs_fo, fam),width = 14, height = 6)
      for (ogi in list_ogs) {
        
        # count architectures per species
        ari = ard[ ard$orthogroup == ogi , ]
        arm = data.frame(xtabs(formula = ~ simarq + species, data=ari))
        arm = tidyr::spread(arm, species, Freq)
        rownames(arm) = arm[,1]
        arm = arm [,-1]
        # order architectures
        arm_gene_counts = colSums(arm)
        arm = rbind(arm_gene_counts, arm)
        rownames(arm)[1] = "Gene presence"
        arm_rowlabs = stringr::str_trunc(rownames(arm), width = 60, ellipsis = "...")
        arm_rowlabs = stringr::str_pad(arm_rowlabs, width = 60, side = "right")
        arm_ord = order(rowSums(arm>0) , apply(arm, 1, which.max), decreasing = T)
        
        # plot duets heatmap
        arp_wte_bool = grepl(list_te_doms_rexp, rownames(arm[arm_ord,]))
        pheatmap(arm[arm_ord,], 
                 color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
                 cellwidth = 4, cellheight = 4, na_col = "grey", gaps_row = 1,
                 cluster_rows = F, cluster_cols = F,
                 fontsize = 5,
                 annotation_row = data.frame(is_TE = arp_wte_bool*1, row.names = rownames(arm[arm_ord,])),
                 gaps_col = sps_gap_ixs, 
                 labels_row = arm_rowlabs[arm_ord],
                 main = sprintf("%s duets, OG %s", fam, ogi),
                 border_color = "white", display_numbers = F, number_format = "%i")
        
        # if there is a TE in this OG, store it for a later concatenated table
        if (any(arp_wte_bool)) {
          
          arm_f = arm[arm_ord,]
          # arm_f = arm_f[rownames(arm_f) != "Gene presence",]
          rownames(arm_f) = paste(
            stringr::str_trunc(ogi, width = 40, ellipsis = "..."),  
            rownames(arm_f), sep = " | ")
          arm_cat = rbind(arm_cat, arm_f)
          arm_cat_breaks = c(arm_cat_breaks, nrow(arm_cat))
          
        }
      }
      dev.off()    
      
    }
  } else { 
    print(sprintf("Can't find %s/net.%s.network_genes.csv!",arqs_fo, fam))  
  }
}


#### Gene trees with TEs ####
# plot all gene trees that have TE insertions, and the protein domain architectures 
# of each gene (simplified)

# bottom, left, top, right
# par(mar = c(5.1, 4.1, 4.1, 2.1))
for (fam in list_fams) {
  
  if (file.exists(sprintf("%s/net.%s.network_genes.csv",arqs_fo, fam))) {
    
    # load arqs for that gene family
    arq = read.table(sprintf("%s/net.%s.network_genes.csv",arqs_fo, fam), sep = "\t", stringsAsFactors = F, header = T)
    arq = arq[arq$family == fam,]
    # simplify architectures (collapse **consecutive** repeats)
    simple_arqs_list = stringr::str_split(arq$architecture, pattern = " ")
    simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
    arq$simarq = simple_arqs_vect
    
    arq_w_TE = grep(list_te_doms_rexp, arq$simarq, perl = T)
    arq_f = arq[arq_w_TE,]
    
    if (nrow(arq_f) > 0)  {
      pdf(file = sprintf("%s/tree_with_TE.%s.pdf", arqs_fo, fam),width = 12, height = 12)
      par(mar = c(1, 1, 1, 1))
      
      ogs_w_TE_fam = stringr::str_split(arq_f$orthogroup, pattern = "\\.", simplify = T)[,1]
      ogs_w_TE_hom = stringr::str_split(arq_f$orthogroup, pattern = "\\.", simplify = T)[,2]
      ogs_w_TE = unique(paste(ogs_w_TE_fam,ogs_w_TE_hom, sep="."))
      
      for (ogi in ogs_w_TE) {
        
        if (file.exists(sprintf("%s/euk.%s.seqs.iqtree.treefile", tree_fo, ogi))) {
          
          gphy = read.tree(
            file=sprintf("%s/euk.%s.seqs.iqtree.treefile", tree_fo, ogi))
          
          # shorten tip labels
          gphy$tip.label = stringr::str_split(gphy$tip.label, pattern = "\\|", simplify = T)[,1]
          
          # which tips have TEs?    
          arq_fi = arq[arq$gene %in% gphy$tip.label, ]
          rownames(arq_fi) = arq_fi$gene
          arq_fi = arq_fi[gphy$tip.label,]
          tips_w_TE = grep(list_te_doms_rexp, arq_fi$simarq)
          
          # add architectures to label
          arqlabels = paste(gphy$tip.label, arq_fi$simarq, sep ="  |  ")
          gphy$tip.label = arqlabels
          
          # plot gene tree      
          text_size=0.3
          plot.phylo(gphy, font=1, type="u", label.offset = 0.05, cex=text_size, edge.color = "gray",
                     root.edge = T, show.tip.label = T, align.tip.label = F, underscore = T,tip.color = alpha("gray10",0.7),
                     main = sprintf("%s\nn = %i genes with TE domains", ogi, length(tips_w_TE)))
          # plot presences
          tiplabels(pch = 19, col = "blue", height = 4, cex = 0.5, tip = tips_w_TE)
          nodelabels(gphy$node.label, col = alpha("darkgray",0.6), frame = "none", cex=text_size)
          # add scale bar
          add.scale.bar(x=0, y=-5)
          
        }
        
      }
      dev.off()
    }
    
  } else {
    print(sprintf("Can't find %s/net.%s.network_genes.csv!",arqs_fo, fam))  
  }
}



### global architecture analyses ###
garq = read.table(garq_fn, header = F, sep="\t", col.names = c("gene","architecture"))
gcla = read.table(gcla_fn, header = F, sep = "\t", col.names = c("gene","gene_fam"))
garq_cla = merge(garq, gcla, by.x = "gene", by.y="gene", all.x = F, all.y = T)

# split architectures into per-domain entries
qarq_list = stringr::str_split(garq_cla$architecture, pattern = " ")
names(qarq_list) = paste(garq_cla$gene, "num")
qarq_lidf = data.frame(domain = unlist(qarq_list, use.names = F))
qarq_lidf$gene = str_split(names(unlist(qarq_list)), " ", simplify = T)[,1]

# reduce DF to TE list
ixs_lidf_wTE = grep(pattern = list_te_doms_rexp, qarq_lidf$domain)
qarq_lidf_wTE = qarq_lidf[ixs_lidf_wTE,]
# qarq_lidf_wTE$domain = droplevels(qarq_lidf_wTE$domain)

# add classification
qarq_lidf_wTE_cla = merge(qarq_lidf_wTE, gcla, by.x = "gene", by.y="gene", all.x = T, all.y = T)
arm = data.frame(xtabs(formula = ~ domain + gene_fam, data=qarq_lidf_wTE_cla))
arm = tidyr::spread(arm, gene_fam, Freq)
rownames(arm) = arm[,1]
arm = arm [,-1]

# number of TE insertions per family
pdf(file = sprintf("%s/genes_with_TE.barplots.pdf", arqs_fo),height = 8, width = 4)
par(mar = c(5.1, 6.1, 4.1, 2.1))
barplot(as.matrix(arm), horiz = T, cex.names = 0.5, cex.axis = 0.5,xlab = "Frequency",
        col = rainbow(nrow(arm), v = 0.8), las =1, border = NA, main="TE fusions per family")
for (fam in colnames(arm)) {
  arm_fam = arm[,fam]
  names(arm_fam) = paste(rownames(arm), "n =", arm_fam)
  barplot(arm_fam,horiz = T, cex.names = 0.5, cex.axis = 0.5,xlab = "Frequency",
          col = rainbow(nrow(arm), v = 0.8), las =1, border = NA, main = fam)
}

dev.off()



#### global distribution of functional classes ####

garq_cla = merge(garq_cla, cla, by = "gene_fam", all.x = TRUE, all.y = FALSE)

pdf(file = sprintf("%s/genes_global_counts.pdf", arqs_fo),height = 6, width = 12)

# general counts
keep_types = c("Acetylase","Deacetylase", "Methylase","Demethylase","Readers","Remodeller","Chaperones")
garq_cla_f = garq_cla [ garq_cla$gene_type %in% keep_types, ]
par(mar = c(12, 6.1, 4.1, 2.1))
garq_cla_t = table(garq_cla_f$gene_type)[keep_types]
barplot(garq_cla_t, names.arg = sprintf("%s n=%i", names(garq_cla_t),garq_cla_t),
        las = 2, main = sprintf("simple counts N=%i", nrow(garq_cla_f)))

# combinations of classes
count_cla_t = aggregate(gene_type ~ gene, garq_cla_f, function(i) paste(unique(sort(i)), collapse = ","))
count_cla_c = sort(table(count_cla_t$gene_type))
barplot(rev(count_cla_c), names.arg = sprintf("%s n=%i", rev(names(count_cla_c)),rev(count_cla_c)),
        las = 2, main = sprintf("combinations N=%i", nrow(count_cla_t)))

par(mar = c(5.1, 6.1, 4.1, 2.1))
pie(count_cla_c, labels = sprintf("%s n=%i", names(count_cla_c),count_cla_c),
    sub = sprintf("combinations N=%i", nrow(count_cla_t)), col = rainbow(length(count_cla_c), s = 0.9))

# catalytic + reader overlap
garq_cla_f$is_catalytic = garq_cla_f$gene_type %in% c("Acetylase","Deacetylase", "Methylase","Demethylase","Remodeller","Chaperones")
count_cla_t = aggregate(is_catalytic ~ gene, garq_cla_f, function(i) paste(unique(sort(i)), collapse = ","))
count_cla_c = sort(table(count_cla_t[,2]))
par(mar = c(5.1, 6.1, 4.1, 2.1))
pie(count_cla_c, labels = sprintf("%s n=%i", names(count_cla_c),count_cla_c),
    main = "is catalytic?",
    sub = sprintf("combinations N=%i", nrow(count_cla_t)), col = rainbow(length(count_cla_c), s = 0.9))

# combinations of catalytic stuff with readers
garq_cla_f$catalytic_reader_doms = NA
garq_cla_f [  garq_cla_f$is_catalytic, ]$catalytic_reader_doms = "Catalytic"
garq_cla_f [ !garq_cla_f$is_catalytic, ]$catalytic_reader_doms = garq_cla_f [ !garq_cla_f$is_catalytic, ]$gene_fam
count_cla_t = aggregate(catalytic_reader_doms ~ gene, garq_cla_f, function(i) paste(unique(sort(i)), collapse = ","))
count_cla_c = sort(table(count_cla_t[,2]))
count_cla_c = count_cla_c [ grepl(",",names(count_cla_c)) & grepl("Catalytic",names(count_cla_c)) ]
par(mar = c(12, 6.1, 4.1, 2.1))
barplot(rev(count_cla_c), names.arg = sprintf("%s n=%i", rev(names(count_cla_c)),rev(count_cla_c)),
        las = 2, main = sprintf("overlaps reader + catalytic N=%i", sum(count_cla_c)))


par(mar = c(5.1, 6.1, 4.1, 2.1))
pie(count_cla_c, labels = sprintf("%s n=%i", names(count_cla_c),count_cla_c),
    main = "overlaps reader + catalytic",
    sub = sprintf("combinations N=%i", sum(count_cla_c)), col = rainbow(length(count_cla_c), s = 0.9))



# combinations of readers
garq_cla_f = garq_cla [ garq_cla$gene_type %in% "Readers", ]
count_cla_t = aggregate(gene_fam ~ gene, garq_cla_f, function(i) paste(unique(sort(i)), collapse = ","))
count_cla_c = sort(table(count_cla_t$gene_fam))
par(mar = c(12, 6.1, 4.1, 2.1))
barplot(rev(count_cla_c), names.arg = sprintf("%s|n=%i", rev(names(count_cla_c)),rev(count_cla_c)),
        las = 2, main = sprintf("combinations readers N=%i", nrow(count_cla_t)))

barplot(rev(count_cla_c[ grep(",",names(count_cla_c))]), names.arg = sprintf("%s|n=%i", rev(names(count_cla_c[ grep(",",names(count_cla_c))])),rev(count_cla_c[ grep(",",names(count_cla_c))])),
        las = 2, main = sprintf("combinations readers N=%i", nrow(count_cla_t)))

par(mar = c(5.1, 6.1, 4.1, 2.1))
pie(count_cla_c, labels = sprintf("%s n=%i", names(count_cla_c),count_cla_c),
    sub = sprintf("combinations readers N=%i", nrow(count_cla_t)), col = rainbow(length(count_cla_c), s = 0.9))

# sum which readers have combinations
ix_combinations = grepl(",",names(count_cla_c))
nu_combinations = sum(count_cla_c[ix_combinations])
ix_pure = !grepl(",",names(count_cla_c))
nu_pure = sum(count_cla_c[ix_pure])
pie(c(nu_pure, nu_combinations), labels = c(paste("pure", nu_pure), paste("combinations", nu_combinations)), main = "readers combinations")



dev.off()
