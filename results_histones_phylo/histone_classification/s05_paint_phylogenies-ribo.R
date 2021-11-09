# libraries
library(stringr)
library(ape)
library(phangorn)
library(data.table)
library(scales)
library(adephylo)

#### concatenated tree ####

# this tree of a concatenation of CCs manually selected because they include canonical 
# eukaryotic histones, as well as viral histone homologs

# input
tre_fn = "alignments/ribo.concatenated/cat_ribohist.iqtree.alrt.treefile"
# tables: virus taxonomy
vtax_fn = "../../data/seq_Viruses.taxa.csv.gz"    # taxon for each viral seq
taxne_fn = "../../data/taxonomy_ranked_wg.tsv.gz" # taxonomy for non-euks
# histone classification from blast
hcl_fn = "all.Histone.to_histdb.Spring-Naive.csv"
# list of outgroup sequences (archaea)
oar_fn = "alignments/Erives_2017.archaeal_histones.to_NCBI.diamond.hits_nr.c95.fasta"

# load outgroup sequences
oar = read.FASTA(oar_fn)
outgroup_archaea = names(oar)
outgroup_archaea = unique(gsub(" .*", "", outgroup_archaea))

# load tree
tre = ape::read.tree(tre_fn)
# root
# tre = phangorn::midpoint(tre)
tre = ape::root(phy = tre, outgroup=outgroup_archaea)
tre = ape::ladderize(tre)

# load histone classification
hcl = read.table(hcl_fn, sep="\t", header = T)
hcl = hcl[!grepl("^vir_|^arc_",hcl$members),]

# find tips
tips_H2A_ix = which(tre$tip.label %in% hcl[hcl$classification=="H2A","members"])
tips_H2B_ix = which(tre$tip.label %in% hcl[hcl$classification=="H2B","members"])
tips_H3_ix = which(tre$tip.label %in% hcl[hcl$classification=="H3","members"])
tips_H4_ix = which(tre$tip.label %in% hcl[hcl$classification=="H4","members"])
tips_H2AZ_ix = which(tre$tip.label %in% hcl[hcl$classification=="H2AZ","members"])
tips_H2Am_ix = which(tre$tip.label %in% hcl[hcl$classification=="macroH2A","members"])
tips_H3CE_ix = which(tre$tip.label %in% hcl[hcl$classification=="other_CENP","members"])
cahis_ixs = c(tips_H3_ix, tips_H4_ix, tips_H2A_ix, tips_H2B_ix, tips_H2AZ_ix, tips_H2Am_ix, tips_H3CE_ix)
other_ixs = which(!(1:length(tre$tip.label) %in% cahis_ixs))

# archaeal sequences
tips_arc_ix = which(grepl("^arc_", tre$tip.label))
tips_vir_ix = which(grepl("^vir_", tre$tip.label))
tips_hsa_ix = which(grepl("^Hsap_", tre$tip.label))

# read data
# taxonomy of viral sequences
vtax = fread(cmd=sprintf("zcat %s", vtax_fn), data.table = F, header = F, stringsAsFactors = F)
colnames(vtax) = c("gene", "taxon")
vtax$gene = paste("vir_", vtax$gene, sep="")
# taxonomy of species
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=F, stringsAsFactors=F)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")

# add taxonomic info to tips (only for viruses)
phyl_nodt = data.frame(gene = tre$tip.label, gene_id = gsub("_\\d+-\\d+$", "", tre$tip.label))
phyl_nodt$gene_id = gsub("_\\d+-\\d+\\|Erives2017$", "", phyl_nodt$gene_id)
phyl_nodt = merge(phyl_nodt, vtax, by.x = "gene_id", by.y = "gene", all.x = T, all.y=F)
phyl_nodt = merge(phyl_nodt, taxne, by.x = "taxon", by.y = "taxon", all.x = T, all.y = F)
# reorder
phyl_nodt = phyl_nodt[match(tre$tip.label, phyl_nodt$gene),]

# pairwise distances
#phyl_dist = adephylo::distTips(tre, method = "patristic")
#phyl_dism = as.matrix(phyl_dist)
#ix_vir = which(grepl("^vir_", rownames(phyl_dist)))
#ix_cel = which(!grepl("^vir_", rownames(phyl_dist)))
#phyl_dism_vir = phyl_dism[ix_vir,ix_cel]
#phyl_dism_vir_closest_ix = apply(phyl_dism_vir, 1, which.min)
#phyl_dism_vir_closest_di = apply(phyl_dism_vir, 1, min)
#phyl_dism_vir_dat = data.frame(
#    gene = rownames(phyl_dism),
#    dist_to_closest = phyl_dism_vir_closest_di,
#    closest_gene = colnames(phyl_dism)[phyl_dism_vir_closest_di]
#)

# plot
pdf(sprintf("../trees_concatenated_ribohistone_CC.reduced.pdf", tre_fn), width = 8, height = 8)
ape::plot.phylo(tre, underscore = T, cex = 0.3, font=1, type = "u", show.tip.label =  F, edge.color = "gray", show.node.label = T)
ape::tiplabels(pch = 16, tip=tips_H3CE_ix, col = "cyan3", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2AZ_ix, col = "chartreuse3", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2Am_ix, col = "darkgreen", cex=0.5)
ape::tiplabels(pch = 16, tip=tips_H3_ix, col = "blue", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H4_ix, col = "orange", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2A_ix, col = "springgreen4", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2B_ix, col = "purple", cex = 0.5)
ape::tiplabels(pch = 19, tip=other_ixs, col = "darkgray", cex = 0.5)
ape::tiplabels(pch = 19, col = "red", cex = 0.5, tip = tips_arc_ix)
ape::tiplabels(pch = 19, col = "magenta", cex = 0.3, tip = tips_vir_ix)
ape::tiplabels(
    text = paste(phyl_nodt$family, phyl_nodt$gene_id)[tips_vir_ix],
    col = alpha("darkmagenta",0.6), 
    tip = tips_vir_ix,
    frame="none", cex=0.4)
ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
legend(
    "topright", 
    legend = c("H3","H4","H2A","H2B","H2AZ","macroH2A","Archaea","Virus"),
    col = c("blue","orange","springgreen4","purple","chartreuse3","darkgreen","red","magenta"), pch=19, cex=0.5)
dev.off()

# plot
pdf(sprintf("../trees_concatenated_ribohistone_CC.pdf", tre_fn), width = 40, height = 40)
ape::plot.phylo(tre, underscore = T, cex = 0.3, font=1, type = "u", show.tip.label =  F, edge.color = "gray", show.node.label = T)
ape::tiplabels(pch = 16, tip=tips_H3CE_ix, col = "cyan3", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2AZ_ix, col = "chartreuse3", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2Am_ix, col = "darkgreen", cex=0.5)
ape::tiplabels(pch = 16, tip=tips_H3_ix, col = "blue", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H4_ix, col = "orange", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2A_ix, col = "springgreen4", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2B_ix, col = "purple", cex = 0.5)
ape::tiplabels(pch = 19, tip=other_ixs, col = "darkgray", cex = 0.5)
ape::tiplabels(pch = 19, col = "red", cex = 0.5, tip = tips_arc_ix)
ape::tiplabels(pch = 19, col = "magenta", cex = 0.3, tip = tips_vir_ix)
ape::tiplabels(
    text = paste(phyl_nodt$family, phyl_nodt$gene_id)[tips_vir_ix],
    col = alpha("darkmagenta",0.6), 
    tip = tips_vir_ix,
    frame="none", cex=0.4)
ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
legend(
    "topright", 
    legend = c("H3","H4","H2A","H2B","H2AZ","macroH2A","Archaea","Virus"),
    col = c("blue","orange","springgreen4","purple","chartreuse3","darkgreen","red","magenta"), pch=19, cex=0.5)
dev.off()


# plot
pdf(sprintf("../trees_concatenated_ribohistone_CCp.pdf", tre_fn), width = 40, height = 200)
ape::plot.phylo(tre, underscore = T, cex = 0.3, font=1, type = "p", show.tip.label =  T, edge.color = "gray", show.node.label = T)
ape::tiplabels(pch = 16, tip=tips_H3CE_ix, col = "cyan3", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2AZ_ix, col = "chartreuse3", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2Am_ix, col = "darkgreen", cex=0.5)
ape::tiplabels(pch = 16, tip=tips_H3_ix, col = "blue", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H4_ix, col = "orange", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2A_ix, col = "springgreen4", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2B_ix, col = "purple", cex = 0.5)
ape::tiplabels(pch = 19, tip=other_ixs, col = "darkgray", cex = 0.5)
ape::tiplabels(pch = 19, col = "red", cex = 0.5, tip = tips_arc_ix)
ape::tiplabels(pch = 19, col = "magenta", cex = 0.3, tip = tips_vir_ix)
ape::tiplabels(
    text = paste(phyl_nodt$family, phyl_nodt$gene_id)[tips_vir_ix],
    col = alpha("darkmagenta",0.6), 
    tip = tips_vir_ix,
    frame="none", cex=0.4)
ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
legend(
    "topright", 
    legend = c("H3","H4","H2A","H2B","H2AZ","macroH2A","Archaea","Virus"),
    col = c("blue","orange","springgreen4","purple","chartreuse3","darkgreen","red","magenta"), pch=19, cex=0.5)
dev.off()


