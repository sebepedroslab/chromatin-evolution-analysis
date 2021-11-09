# libraries
library("stringr")
library("ape")
library("phangorn")
library("data.table")
library("scales")
library("phytools")

# input: eukaryotic-only phylogenies
tre_li = c("alignments/euk.Histone.to_histdb/CC4_H2A.iqt.treefile",
           "alignments/euk.Histone.to_histdb/CC0_H2B.iqt.treefile",
           "alignments/euk.Histone.to_histdb/CC2_H3.iqt.treefile",
           "alignments/euk.Histone.to_histdb/CC1_H4.iqt.treefile")
his_li = c("H2A","H2B","H3","H4")
col_li = c("springgreen3","purple","blue","orange")
out_fo = "histone_trees/"

# histones that appear in proteomics
pro_fn = "all.Histone.to_proteomics.csv.gz" # histones that appear in the proteomics datasets
mod_fn = "../../results_PTMs/consensus_modifications_perseq.xlsx" # list of canonical histones used in proteomics analyse

# load list of histones that appear in the proteomics dataset
pro = read.table(pro_fn, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
# keep only identical hits
pro = pro [ pro$pident == 100, ]

# list of canonical histones identified in the proteomics dataset
# keep only reference sequences that appear in the proteomics dataset
his_in_pro = list()
hip_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")
hip_cols = c("blue4","darkorange3","darkgreen","darkmagenta","darkolivegreen3","darkolivegreen4")
for (hip in hip_list) {
     
     mod = read_excel(mod_fn, sheet=hip)
     mod = mod [ !is.na(mod$Positions_in_Reference), ]
     his_in_mod = gsub("\\|.*","",mod$Protein_ID)
     his_in_pro[[hip]] = pro [ pro$sseqid %in% his_in_mod , "qseqid" ]
     
}

pdf("../trees_histone_CC.pdf", height = 10, width = 5)
layout(matrix(1:4, ncol = 2))
for (n in seq_along(tre_li)) {
    
    tre_fn = tre_li[n] 
    tre = ape::read.tree(tre_fn)
    tre = phangorn::midpoint(tre)    
    tre = ape::ladderize(tre) 
    ape::plot.phylo(tre, edge.color = col_li[n], underscore = TRUE, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  FALSE, main=his_li[n])
    ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
    
    for (m in 1:length(hip_list)) {
        hip = hip_list[m]
        hic = hip_cols[m]
        ix_hip = which(tre$tip.label %in% his_in_pro[[hip]])
        if (length(ix_hip) > 0) {
            ape::tiplabels(tip = ix_hip, pch=1, col=hic, cex=0.8)
        }
    }
}
dev.off()

pdf("../trees_histone_CC_names.pdf", height = 30, width = 5)
for (n in seq_along(tre_li)) {
    
    tre_fn = tre_li[n] 
    
    tre = ape::read.tree(tre_fn)
    tre = phangorn::midpoint(tre)    
    tre = ape::ladderize(tre) 
    ape::plot.phylo(tre, edge.color = col_li[n], underscore = TRUE, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  TRUE, main=his_li[n])
    ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
}
dev.off()


# input: all histones (includes archaea and viruses)
tre_li = c("alignments/all.Histone.to_histdb/CC34_H2A.iqt.treefile",
           "alignments/all.Histone.to_histdb/CC33_H2B.iqt.treefile",
           "alignments/all.Histone.to_histdb/CC30_H3.iqt.treefile",
           "alignments/all.Histone.to_histdb/CC31_H4.iqt.treefile")
his_li = c("H2A","H2B","H3","H4")
col_li = c("springgreen3","purple","blue","orange")
out_fo = "histone_trees/"
seq_li = list()

pdf("../trees_histone_CCnoneuk.pdf", height = 10, width = 5)
layout(matrix(1:4, ncol = 2))
for (n in seq_along(tre_li)) {
    
    tre_fn = tre_li[n] 
    
    tre = ape::read.tree(tre_fn)
    tre = phangorn::midpoint(tre)
    tre = ape::ladderize(tre)
    ix_vir = which(grepl("^vir_",tre$tip.label ))
    ix_arc = which(grepl("^arc_",tre$tip.label ))
    ape::plot.phylo(tre, edge.color = col_li[n], underscore = TRUE, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  FALSE,
                    main=his_li[n])
    ape::tiplabels(tip = ix_vir, pch=1, col="red", cex=0.8)
    ape::tiplabels(tip = ix_arc, pch=1, col="darkorange", cex=0.8)
    ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
    
    # store sequence list
    seq_li[[n]] = tre$tip.label
    
}
dev.off()

pdf("../trees_histone_CCnoneuk_names.pdf", height = 30, width = 5)
for (n in seq_along(tre_li)) {
    
    tre_fn = tre_li[n] 
    
    tre = ape::read.tree(tre_fn)
    tre = phangorn::midpoint(tre)
    tre = ape::ladderize(tre)
    ape::plot.phylo(tre, edge.color = col_li[n], underscore = TRUE, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  TRUE,
                    main=his_li[n])
    ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)

}

dev.off()



#### concatenated tree ####

# this tree of a concatenation of CCs manually selected because they include canonical 
# eukaryotic histones, as well as viral histone homologs

# input
tre_fn = "alignments/vir.concatenated/cat_histones.iqtree.alrt.treefile"
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
hcl = read.table(hcl_fn, sep="\t", header = TRUE)
hcl = hcl[!grepl("^vir_|^arc_",hcl$members),]

# find tips
tips_H2A_ix = which(tre$tip.label %in% hcl[hcl$classification == "H2A","members"])
tips_H2B_ix = which(tre$tip.label %in% hcl[hcl$classification == "H2B","members"])
tips_H3_ix = which(tre$tip.label %in% hcl[hcl$classification == "H3","members"])
tips_H4_ix = which(tre$tip.label %in% hcl[hcl$classification == "H4","members"])
tips_H2AZ_ix = which(tre$tip.label %in% hcl[hcl$classification == "H2AZ","members"])
tips_H2Am_ix = which(tre$tip.label %in% hcl[hcl$classification == "macroH2A","members"])
tips_H3CE_ix = which(tre$tip.label %in% hcl[hcl$classification == "other_CENP","members"])
cahis_ixs = c(tips_H3_ix, tips_H4_ix, tips_H2A_ix, tips_H2B_ix, tips_H2AZ_ix, tips_H2Am_ix, tips_H3CE_ix)
other_ixs = which(!(1:length(tre$tip.label) %in% cahis_ixs))

# archaeal sequences
tips_arc_ix = which(grepl("^arc_", tre$tip.label))
tips_vir_ix = which(grepl("^vir_", tre$tip.label))
tips_hsa_ix = which(grepl("^Hsap_", tre$tip.label))

# read data
# taxonomy of viral sequences
vtax = fread(cmd=sprintf("zcat %s", vtax_fn), data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
colnames(vtax) = c("gene", "taxon")
vtax$gene = paste("vir_", vtax$gene, sep="")
# taxonomy of species
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=FALSE, stringsAsFactors=FALSE)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")

# add taxonomic info to tips (only for viruses)
phyl_nodt = data.frame(gene = tre$tip.label, gene_id = gsub("_\\d+-\\d+$", "", tre$tip.label))
phyl_nodt$gene_id = gsub("_\\d+-\\d+\\|Erives2017$", "", phyl_nodt$gene_id)
phyl_nodt = merge(phyl_nodt, vtax, by.x = "gene_id", by.y = "gene", all.x = TRUE, all.y=FALSE)
phyl_nodt = merge(phyl_nodt, taxne, by.x = "taxon", by.y = "taxon", all.x = TRUE, all.y = FALSE)
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
pdf(sprintf("../trees_concatenated_histone_CC.reduced.pdf", tre_fn), width = 8, height = 8)
ape::plot.phylo(tre, underscore = TRUE, cex = 0.3, font=1, type = "u", show.tip.label =  FALSE, edge.color = "gray", show.node.label = TRUE)
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
pdf(sprintf("../trees_concatenated_histone_CC.pdf", tre_fn), width = 40, height = 40)
ape::plot.phylo(tre, underscore = TRUE, cex = 0.3, font=1, type = "u", show.tip.label =  FALSE, edge.color = "gray", show.node.label = TRUE)
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


#### summary viral histones ####

# from concatenated tree
# create table
virh = phyl_nodt[grepl("^vir", phyl_nodt$gene),c("gene","gene_id","taxon","family","superkingdom")]

# add similarity info:
alis = read.table(
    "all.Histone.to_histdb.csv.gz", header = FALSE, sep = "\t",
    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
    stringsAsFactors = FALSE)
# find a reference canonical histone and record distance from virus to it:
for (hii in c("H2A","H2B","H3","H4")) {
    print(sprintf("dist to %s...", hii))
    ref_his = alis[alis$pident  ==  100 & alis$length > 70 & grepl(sprintf("\\|%s$", hii), alis$sseqid),"qseqid"][1]
    virh[,sprintf("dist_to_%s",hii)] = unlist(lapply(virh$gene, FUN = function(x) phytools::fastDist(tre, x, ref_his)))
}
# which is closest?
closest_canonical_ix = apply(virh[c("dist_to_H2A","dist_to_H2B","dist_to_H3","dist_to_H4")], 1, which.min)
virh$closest_canonical = c("H2A","H2B","H3","H4")[closest_canonical_ix]

# sort by taxon
virh = virh[order(virh$family,virh$taxon,virh$gene_id,virh$closest_canonical),]

# save
write.table(virh, file = "../viral_histones_summary.csv", row.names = FALSE, sep="\t", quote = FALSE)

