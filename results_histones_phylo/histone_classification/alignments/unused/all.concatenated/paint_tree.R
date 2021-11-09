library(ape)
library(phangorn)

# input
tre_fn = "allhistone_dom.c99.ali_untrimmed.iqtf.treefile"
# tre_fn = "allhistone_dom.c99.ali.iqtfb_untrimmed.treefile"

# load tree
tre = ape::read.tree(tre_fn)

# find tips
tips_H3_ix = which(grepl("_H3$", tre$tip.label))
tips_H4_ix = which(grepl("_H4$", tre$tip.label))
tips_H2A_ix = which(grepl("_H2A$", tre$tip.label))
tips_H2B_ix = which(grepl("_H2B$", tre$tip.label))
cahis_ixs = c(tips_H3_ix, tips_H4_ix, tips_H2A_ix, tips_H2B_ix)
other_ixs = which(!(1:length(tre$tip.label) %in% cahis_ixs))
    
tips_arc_ix = which(grepl("^arc_", tre$tip.label))
tips_vir_ix = which(grepl("^vir_", tre$tip.label))
tips_hsa_ix = which(grepl("^Hsap_", tre$tip.label))

# root
tre = phangorn::midpoint(tre)
tre = ape::ladderize(tre)

# plot
pdf(sprintf("%s.pdf", tre_fn), width = 10, height = 500)
ape::plot.phylo(tre, underscore = T, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  F)
ape::tiplabels(text=tre$tip.label[tips_H3_ix], tip=tips_H3_ix, col = "blue", cex = 0.5, bg = F, frame = "none", adj=0 )
ape::tiplabels(text=tre$tip.label[tips_H4_ix], tip=tips_H4_ix, col = "orange", cex = 0.5, bg = F, frame = "none", adj=0 )
ape::tiplabels(text=tre$tip.label[tips_H2A_ix], tip=tips_H2A_ix, col = "springgreen4", cex = 0.5, bg = F, frame = "none", adj=0 )
ape::tiplabels(text=tre$tip.label[tips_H2B_ix], tip=tips_H2B_ix, col = "purple", cex = 0.5, bg = F, frame = "none", adj=0 )
ape::tiplabels(text=tre$tip.label[other_ixs], tip=other_ixs, col = "darkgray", cex = 0.5, bg = F, frame = "none", adj=0 )
ape::tiplabels(pch = 19, col = "magenta", height = 4, cex = 0.5, tip = tips_vir_ix, offset=-0.02)
ape::tiplabels(pch = 19, col = "red", height = 4, cex = 0.5, tip = tips_hsa_ix, offset=-0.02)
dev.off()



# plot
pdf(sprintf("%s.reduced.pdf", tre_fn), width = 20, height = 20)
ape::plot.phylo(tre, underscore = T, cex = 0.2, font=1, col="grey",type = "u", show.tip.label =  F)
ape::tiplabels(pch = 19, tip=tips_H3_ix, col = "blue", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H4_ix, col = "orange", cex = 0.5)
ape::tiplabels(pch = 19, tip=tips_H2A_ix, col = "springgreen4", cex=0.5)
ape::tiplabels(pch = 19, tip=tips_H2B_ix, col = "purple", cex = 0.5)
ape::tiplabels(pch = 19, tip=other_ixs, col = "darkgray", cex = 0.5)
ape::tiplabels(pch = 19, col = "magenta", cex = 0.5, tip = tips_vir_ix)
dev.off()
