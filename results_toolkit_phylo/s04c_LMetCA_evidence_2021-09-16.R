#### Input ####
# load libraries
library(ape)
library(stringr)

pres_fn = "orthogroups_euk.ancestral.posteriors_pres.csv"
gain_fn = "orthogroups_euk.ancestral.posteriors_gain.csv"
summ_fn = "orthogroups_euk.ancestral.posteriors_summ.csv"
dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
outp_fn = "results_evol/"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"

prob_thr = 0.5

graphics.off()

#### Define input ####

# load data
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE)
cla = read.table(clas_fn, header = FALSE, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_counts = table(dat$gene)
genes_repeated = names(gene_counts[gene_counts > 1])
ogs_lis_w_dup_genes = dat[dat$gene %in% genes_repeated,"orthogroup"]
dat_dups_full = dat[dat$orthogroup %in% ogs_lis_w_dup_genes,]

# list ogs
ogs_lis = unique(dat_dups_full$orthogroup)
ogs_lis_black = c()

# read species tree
phyl = ape::read.tree(file = phyl_fn)
phyl$edge.length = rep(1, nrow(phyl$edge))
# phyl = phytools::rotateNodes(phyl, nodes = "all")
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
phyl_edge = merge(phyl_edge, phyl_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")


# summary of OG data
mat_summ = read.table(summ_fn, sep="\t", header = TRUE)
# remove blacklisted OGs (redundant with others)
mat_summ = mat_summ[!(mat_summ$orthogroup %in% ogs_lis_black),]

# fix OG ids and classes
mat_summ$og_class = stringr::str_split(mat_summ$orthogroup,pattern = "\\.", simplify = TRUE)[,1]
mat_summ$og_id    = stringr::str_split(mat_summ$orthogroup,pattern = ":", simplify = TRUE)[,1]
# fix OG names
mat_summ$og_tmp   = stringr::str_remove(mat_summ$orthogroup, pattern = ":likeclu:.*")
mat_summ$og_tmp   = stringr::str_replace(mat_summ$og_tmp,    pattern = ":like:", replacement = ":like_")
mat_summ$og_name  = stringr::str_split(mat_summ$og_tmp,pattern = ":", simplify = TRUE)[,2]
mat_summ$og_name  = stringr::str_replace(mat_summ$og_name,    pattern = "^like_", replacement = "like:")
mat_summ$og_name_short = stringr::str_trunc(mat_summ$og_name, width = 50)
mat_summ$og_id_short   = paste(mat_summ$og_id, mat_summ$og_name_short, sep=" | ")
rownames(mat_summ) = mat_summ$og_id
mat_summ_original_order = rownames(mat_summ)
# add enzyme classes
# mat_summ = merge(mat_summ, cla, by.x = "og_class", by.y = "gene_fam", all.x = TRUE, all.y = FALSE)
# rownames(mat_summ) = mat_summ$og_id
# mat_summ = mat_summ[mat_summ_original_order,]

# load ancestral presences
mat_pres = read.table(pres_fn, sep="\t", header = TRUE, row.names = 1)
mat_gain = read.table(gain_fn, sep="\t", header = TRUE, row.names = 1)

# add num species presence (binarised)
mat_summ$n_presence_bin = apply(mat_pres[,phyl$tip.label], MARGIN = 1, function(x) sum(x > 0))

# simplify OG names
rownames(mat_pres) = mat_summ$og_id
rownames(mat_gain) = mat_summ$og_id

# remove blacklisted OGs
print(sprintf("Remove %i OGs", length(ogs_lis_black)))
print(ogs_lis_black)
mat_summ = mat_summ[!(rownames(mat_summ) %in% ogs_lis_black),]
mat_pres = mat_pres[!(rownames(mat_pres) %in% ogs_lis_black),]
mat_gain = mat_gain[!(rownames(mat_gain) %in% ogs_lis_black),]

# order taxa according to phylogeny?
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]

# create factors for taxa and families
# mat_summ$gain = factor(mat_summ$gain, levels = taxa_lis)
mat_summ$og_class = factor(mat_summ$og_class, levels = unique(mat_summ$og_class))

# add classes to per-gene table
dat = merge(dat, mat_summ, by.x="orthogroup", by.y="orthogroup", all.x = TRUE, all.y = FALSE)
# add species 
dat$taxa = stringr::str_split(dat$gene,pattern = "_", simplify = TRUE)[,1]
dat$taxa = factor(dat$taxa, levels = taxa_lis)
# remove blacklisted ogs
dat = dat[!(dat$og_class %in% ogs_lis_black), ]

# define colors heatmap
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","dodgerblue4"))

# list of gene types
gene_types = unique(as.character(mat_summ$gene_type))



#### Reconstruct LMCA ####

# restrict analysis to OGs present in LMCA
mat_gain [ is.na(mat_gain) ] = 0
mat_summ_i = mat_summ[mat_gain$Metazoa > prob_thr ,]
mat_pres_i = mat_pres[mat_gain$Metazoa > prob_thr ,]

# binarise presence heatmap
# mat_pres_i = (mat_pres_i > 0) * 1

# restrict presence table to taxa that directly descend from LMCA
leca_descendants = as.character(phyl_edge[phyl_edge$edge_start == length(phyl$tip.label) + which(phyl$node.label == "Metazoa"), "taxa"])
leca_descendants = leca_descendants[!is.na(leca_descendants)]
mat_pres_i = mat_pres_i[,leca_descendants]

# which is the support of each lineage?
mat_pres_i_taxid = apply(mat_pres_i, 1, function(x) which(x > prob_thr))
mat_pres_i_taxid_string = unlist(lapply(mat_pres_i_taxid, function(x) paste(sort(names(x)), collapse = ",")))

# store results
mat_summ_i$ancestral_support_string = mat_pres_i_taxid_string
mat_summ_i$ancestral_support_count = stringr::str_count(mat_pres_i_taxid_string, ",") + 1

# count how many LMCA ogs are supported by which ancestral lineages
pdf("results_evol/lmca_gainreconstruction_evidence.pdf", width = 6, height = 12)
support_string = sort(table(mat_summ_i$ancestral_support_string))
par(mar=c(3.1, 14.1, 4.1, 2.1))
barplot(support_string, las=1, horiz = TRUE, cex.names = 0.5,
        main = "LMCA presence")
dev.off()


# count how many LMCA ogs are supported by good ancestral lineages, or by a certain amount of lineages (ad-hoc)
mat_summ_i$ancestral_support_summary = mat_summ_i$ancestral_support_string
mat_summ_i$ancestral_support_summary [ mat_summ_i$ancestral_support_summary == "" ] = "unk"

# factor
mat_summ_i$ancestral_support_summary = factor(mat_summ_i$ancestral_support_summary, levels = unique(mat_summ_i$ancestral_support_summary))



# save
write.table(mat_summ_i, file = "results_evol/lmca_gainreconstruction_evidence.csv", sep="\t", quote = FALSE, row.names = FALSE)

#### Plot pies ####

# open device
pdf("results_evol/lmca_gainreconstruction_pies_all.pdf",width = 6.5, height = 6.5)
layout(matrix(1:16, nrow = 4))
# general plot
support_summary = table(mat_summ_i$ancestral_support_summary)
pie(support_summary, las=1, 
    col = rainbow(nlevels(mat_summ_i$ancestral_support_summary)),
    cex=0.6, main="all", cex.main=0.8, cex.sub=0.8,
    labels = sprintf("%s\nn=%i", names(support_summary),support_summary),
    sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
pie(1, col = NA, border = NA, labels=NA)
legend("topright", legend=levels(mat_summ_i$ancestral_support_summary), fill=rainbow(nlevels(mat_summ_i$ancestral_support_summary)), cex=0.6)
dev.off()

# now per-class or per-type piecharts

# list of gene types
gene_types = unique(as.character(mat_summ$gene_type))
gene_classes = unique(as.character(mat_summ$gene_class))

# loop for classes
pdf("results_evol/lmca_gainreconstruction_pies_geneclass.pdf",width = 6.5, height = 6.5)
layout(matrix(1:16, nrow = 4))
for (type in gene_classes) {
  support_summary = table(mat_summ_i[mat_summ_i$gene_class == type,]$ancestral_support_summary)
  if (sum(support_summary) > 0) {
    pie(support_summary, las=1, 
        col = rainbow(nlevels(mat_summ_i$ancestral_support_summary)),
        cex=0.6, main=sprintf("Class: %s",type), cex.main=0.8,cex.sub=0.8,
        labels = sprintf("%s\nn=%i", names(support_summary),support_summary),
        sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
  } else {
    plot(NA,NA, main=sprintf("Class: %s\nnot in LMCA", type), axes = FALSE, xlab=NA, ylab=NA, xlim = c(0,0), ylim=c(0,0))
  }
}
pie(1, col = NA, border = NA, labels=NA)
legend("topright", legend=levels(mat_summ_i$ancestral_support_summary), fill=rainbow(nlevels(mat_summ_i$ancestral_support_summary)), cex=0.6)
dev.off()

pdf("results_evol/lmca_gainreconstruction_pies_genetype.pdf",width = 6.5, height = 6.5)
layout(matrix(1:16, nrow = 4))
for (type in gene_types) {
  support_summary = table(mat_summ_i[mat_summ_i$gene_type == type,]$ancestral_support_summary)
  if (sum(support_summary) > 0) {
    pie(support_summary, las=1, 
        col = rainbow(nlevels(mat_summ_i$ancestral_support_summary)),
        cex=0.6, main=sprintf("Type: %s",type), cex.main=0.8,cex.sub=0.8,
        labels = sprintf("%s\nn=%i", names(support_summary),support_summary),
        sub=sprintf("N=%i at p>=%.2f", sum(support_summary),prob_thr))
  } else {
    plot(NA,NA, main=sprintf("Type: %s\nnot in LMCA", type), axes = FALSE, xlab=NA, ylab=NA, xlim = c(0,0), ylim=c(0,0))
  }
}
pie(1, col = NA, border = NA, labels=NA)
legend("topright", legend=levels(mat_summ_i$ancestral_support_summary), fill=rainbow(nlevels(mat_summ_i$ancestral_support_summary)), cex=0.6)
dev.off()


##### barplots of LMCA ####

# barplots of LMCA content: 

pdf("results_evol/lmca_gainreconstruction_barplots.pdf",width = 6, height = 3.5)
par(mar=c(9,2.1,2.1,1.1))

# num that pass the threshold
mat_summ_i$og_class_factor = factor(mat_summ_i$og_class, levels = as.character(cla$gene_fam))
xtab = table(mat_summ_i$og_class_factor)
barplot(xtab, las = 2, main=sprintf("# LMCA OGs at prob >=%.2f", prob_thr),
        names.arg = sprintf("%s n=%i", names(xtab), xtab), cex.names = 0.5, cex.axis = 0.5, cex.main=0.5)

# num that pass the threshold, by evidence
mat_summ_i$og_class_factor = factor(mat_summ_i$og_class, levels = as.character(cla$gene_fam))
xtab = table(mat_summ_i$og_class_factor, mat_summ_i$ancestral_support_summary)
xtab = t(as.matrix(xtab))
barplot(xtab, las = 2, main=sprintf("# LMCA OGs at prob >=%.2f", prob_thr),
        col = rainbow(nlevels(mat_summ_i$ancestral_support_summary)),
        cex.names = 0.5, cex.axis = 0.5, cex.main=0.5)
legend("topright", legend=levels(mat_summ_i$ancestral_support_summary), fill=rainbow(nlevels(mat_summ_i$ancestral_support_summary)), cex=0.4, bty = "n")


# num that pass the threshold, by evidence and type
mat_summ_i$og_class_factor_type = factor(mat_summ_i$gene_type, levels = unique(as.character(cla$gene_type)))
xtab = table(mat_summ_i$og_class_factor_type, mat_summ_i$ancestral_support_summary)
xtab = t(as.matrix(xtab))
barplot(xtab, las = 2, main=sprintf("# LMCA OGs at prob >=%.2f", prob_thr),
        col = rainbow(nlevels(mat_summ_i$ancestral_support_summary)),
        cex.names = 0.5, cex.axis = 0.5, cex.main=0.5)
legend("topright", legend=levels(mat_summ_i$ancestral_support_summary), fill=rainbow(nlevels(mat_summ_i$ancestral_support_summary)), cex=0.4, bty = "n")

# sum of probabilities
pres_leca = mat_pres[,"Metazoa"]
pres_leca [ is.na(pres_leca) ] = 0
mat_summ_leca = aggregate(pres_leca, by=list(names = mat_summ$og_class), sum)
rownames(mat_summ_leca) = mat_summ_leca$names
xtab = mat_summ_leca[as.character(cla$gene_fam),"x"]
names(xtab) = as.character(cla$gene_fam)
barplot(xtab, las = 2, main=sprintf("sum probs LMCA OGs", prob_thr),
        names.arg = sprintf("%s n=%.2f", names(xtab), xtab), cex.names = 0.5, cex.axis = 0.5, cex.main=0.5)

dev.off()
