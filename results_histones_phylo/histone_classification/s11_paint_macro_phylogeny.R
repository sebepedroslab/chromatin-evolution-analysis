# libraries
library(stringr)
library(ape)
library(phangorn)
library(data.table)
library(scales)
library(adephylo)

#### concatenated tree ####

# input
tre_fn = "results_macro/mac.Macro.HG1.domains.iqtree.treefile"
ogs_fn = "results_macro/mac.Macro.HG1.domains.iqtree.treefile.ortholog_groups.csv"
vcl_fn = "variants/table_variants.csv"


# read
tre = ape::read.tree(tre_fn)
ogs = read.table(ogs_fn, header = TRUE)
vcl = read.table(vcl_fn, header = TRUE)

# laddeerise
tre = ape::ladderize(tre)

# which tips are macroh2a
tips_macroh2a = which(tre$tip.label %in% ogs$gene [ ogs$orthogroup == "Macro.HG1.3:MACROH2A1/MACROH2A2" ])
tips_histone  = which(tre$tip.label %in% vcl$gene & tre$tip.label %in% ogs$gene [ ogs$orthogroup == "Macro.HG1.3:MACROH2A1/MACROH2A2" ])

# plot
pdf(sprintf("results_macro/tree_macro.pdf", tre_fn), width = 8, height = 8)

# colored
ape::plot.phylo(tre, underscore = T, cex = 0.3, font=1, type = "u", show.tip.label =  FALSE, edge.color = "gray", show.node.label = FALSE, rotate.tree = 90)
ape::tiplabels(pch = 16, tip=tips_macroh2a, col = "chartreuse2", cex = 0.5)
ape::tiplabels(pch = 16, tip=tips_histone,  col = "chartreuse4", cex=0.5)
ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)

# with supports
ape::plot.phylo(tre, underscore = T, cex = 0.3, font=1, type = "u", show.tip.label =  TRUE, edge.color = "gray", show.node.label = TRUE, rotate.tree = 90)
ape::tiplabels(pch = 16, tip=tips_macroh2a, col = "chartreuse2", cex = 0.5)
ape::tiplabels(pch = 16, tip=tips_histone,  col = "chartreuse4", cex=0.5)
ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)

dev.off()


