# libraries
library("stringr")
library("ape")
library("pheatmap")


phyl_fn = "../../results_toolkit_phylo/data/species_tree.newick"
taxe_fn = "../../data/euk_taxonomy_annotated_2020-08-11.csv"   # taxonomic info for eukaryotes
alis_fn = "variants/euk.Histone.domains.to_hmm_database.tab.csv"
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))
arqs_fn = "../../results_toolkit_phylo/results_domains/net.Histone.network_genes.csv"
seqs_fn = "euk.Histone.seqs.fasta"


# sps order from tree
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label

# load domain architectures
arq = read.table(arqs_fn, sep = "\t", header = TRUE)
arq_macro = arq [ grepl("Macro", arq$architecture),  ]

# load sequences
seq = Biostrings::readAAStringSet(seqs_fn)

# match motif for H2Ax in animals
mot_h2ax1 = Biostrings::AAString("SQEY")
mot_h2ax1_m = Biostrings::vmatchPattern(mot_h2ax1, seq)
mot_h2ax1_b = unlist(lapply(1:length(mot_h2ax1_m), function(i) length(mot_h2ax1_m[[i]]) > 0 ) )
# match motif for H2Ax in plants
mot_h2ax2 = Biostrings::AAString("SQEF")
mot_h2ax2_m = Biostrings::vmatchPattern(mot_h2ax2, seq)
mot_h2ax2_b = unlist(lapply(1:length(mot_h2ax2_m), function(i) length(mot_h2ax2_m[[i]]) > 0 ) )
# sequences that match motif
mot_h2ax_v = unique(c(
    names(mot_h2ax1_m) [ mot_h2ax1_b ],
    names(mot_h2ax2_m) [ mot_h2ax2_b ]
))



# load taxonomy and find out where to add spacing in heatmaps
taxe = read.table(taxe_fn, sep = "\t", header = T, stringsAsFactors = T)
taxe = taxe[match(sps_order, taxe$Species),]
sps_gap_ixs = c(1,1+which(diff(as.numeric(taxe$Group))!=0)) - 1

# load hmmscan
dat = read.table(alis_fn, sep = "", col.names = c("target","accession_t","query","accession_q","eval","score","bias","dom_evalue","dom_score","dom_bias","exp","reg","clu","ov","env","dom","rep","inc","description"))
dat = dat [ order(dat$query, dat$score, decreasing = TRUE), ]
# dat$target [ dat$target == "H2A_macroH2Acustom" ] = "H2A_macroH2A"

# first and second hits
first_ixs = which(!duplicated(dat$query))
nonfirst_ixs = which(duplicated(dat$query))
dat_f = dat [ first_ixs, c("query","target","score") ]
dat_s = dat [ nonfirst_ixs, c("query","target","score") ]
dat_s = dat_s [ !duplicated(dat_s$query), ]

# merge first and second
dat_fs = merge(dat_f, dat_s, by = "query", all.x = TRUE, all.y = TRUE)
colnames(dat_fs) = c("query","first_target","first_score","second_target","second_score")

# ratio between first and second hit
dat_fs$score_ratio = dat_fs$first_score / dat_fs$second_score

# add species
dat_fs$species = stringr::str_split(dat_fs$query, "_", simplify = TRUE)[,1]
dat_fs$species = factor(dat_fs$species, levels = sps_order)

# ensure min quality and ratio of identification
dat_fs = dat_fs [ dat_fs$first_score > 70 ,]
dat_fs$first_target [ dat_fs$score_ratio < 1.01 ] = paste(dat_fs$first_target,"putative") [ dat_fs$score_ratio < 1.01 ]

# for each variant, check if distribution of bitscores is extreme, if it is, report as putative
dat_fs$first_score_p = 1
for (var in c("H2A_H2A.W","H2A_H2A.Z","H2A_H2A.X","H2A_macroH2A","H3_H3.3","H3_cenH3")) { 
    var_bool = dat_fs$first_target == var
    dat_fs$first_score_p [ var_bool ] = pnorm(scale(dat_fs$first_score [ var_bool ] ))
    dat_fs$first_target [ dat_fs$first_score_p < 0.05 & var_bool ] = paste(dat_fs$first_target,"putative") [ dat_fs$first_score_p < 0.05 & var_bool ]
}

# ditch canonical
dat_fs$first_target [ grepl("canonical",dat_fs$first_target)  ] = "NA"

# add clean gene name
dat_fs$gene = gsub("_\\d+-\\d+$","",dat_fs$query)


# add h2ax motif presence info
dat_fs$has_H2AX_motif = (grepl("H2A", dat_fs$first_target) | grepl("H2A", dat_fs$second_target)) & dat_fs$gene %in% mot_h2ax_v
dat_fs$first_target [ dat_fs$has_H2AX_motif ] = "H2A_H2A.X"

# add macro domain presence info
dat_fs$has_Macro_motif = dat_fs$first_target != "H2A_macroH2A" & dat_fs$gene %in% arq_macro$gene
dat_fs$first_target [ dat_fs$has_Macro_motif ] = "H2A_macroH2A"

# crosstabulate
dat_f_t = table(dat_fs$first_target, dat_fs$species)
dat_f_t = dat_f_t [order(rownames(dat_f_t), decreasing = TRUE), ]
dat_f_t = dat_f_t [rownames(dat_f_t) != "NA", ]


dat_f_t [ dat_f_t >= 1 ] = 10

pdf(file = sprintf("variants/table_variants.pdf"),width = 14, height = 6)
pheatmap(
    dat_f_t, 
    color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
    cellwidth = 4, cellheight = 4, na_col = "grey", 
    cluster_rows = F, cluster_cols = F,
    fontsize = 5,
    gaps_col = sps_gap_ixs, 
    gaps_row = grep("putative",rownames(dat_f_t),invert = TRUE),
    main = "variants presence",
    border_color = "white", display_numbers = F, number_format = "%i")
dev.off()


# save: discard unclassified
dat_fs_r = dat_fs [ dat_fs$first_target != "NA", ]
# save: add sequences
dat_fs_r$sequence = seq [ dat_fs_r$gene ]
# save: write
write.table(dat_fs_r, "variants/table_variants.csv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
