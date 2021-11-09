# libraries
library(data.table)
library(pheatmap)
library(stringr)
library(ape)

# input files
taxe_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"   # taxonomic info for eukaryotes
gen_list_fn = "../data/gene_families_hmm.csv"               # list of gene familes
search_fn = "gene_counts/"
phylo_fn = "gene_trees/"
outp_fn = "results_counts/"
spsphy_fn = "data/species_tree.newick"

graphics.off()

# heatmap colors
# col.fun = colorRampPalette(interpolate="l",c("azure2","deepskyblue","dodgerblue4"))
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

# load list hmms
gen_list = read.table(gen_list_fn, header = F, stringsAsFactors = T)
colnames(gen_list) = c("Class","Type","Family","Domains","search_thr","inflation","min_phylo_size")
gen_list_gap_ixs = c(1,1+which(diff(as.numeric(gen_list$Type))!=0)) - 1

# function to order matrix according to size and which.max
order_matrix = function(mat, dim=1) {
    
    if (dim == 1) {
        if (length(mat) > 2) {
            arm_ord = order(rowSums(mat>0) , apply(mat, dim, which.max), decreasing = T)
            mao = mat[arm_ord,]
        } else {
            arm_ord = 1:length(mat)
            mao = mat[arm_ord]
        }
		mao_names = rownames(mat)[arm_ord]
		mao = matrix(mao, nrow = nrow(mat))
		rownames(mao) = mao_names
    } else {
        if (length(mat) > 2) {
            arm_ord = order(colSums(mat>0) , apply(mat, dim, which.max), decreasing = T)
            mao = mat[,arm_ord]
        } else {
            arm_ord = 1:length(mat)
            mao = mat[arm_ord]
        }
		mao_names = colnames(mat)[arm_ord]
		mao = matrix(mao, nrow = nrow(mat))
		colnames(mao) = mao_names
    }
	
    return(mao)
}

#### eukaryotes HMM counts ####

# read taxonomy
taxe = read.table(taxe_fn, sep = "\t", header = T, stringsAsFactors = T)

# load sps tree (for row ordering)
sphy = ape::read.tree(spsphy_fn)
taxe = taxe[match(sphy$tip.label, taxe$Species),]

# read gene counts
gen  = fread(input = "gene_counts/euk_genecounts.csv", data.table = F, header = F, col.names = c("gene","gene_family"))
gen$species = stringr::str_split(gen$gene, pattern = "_", simplify = T)[,1]

# gene families as ordered factors
gen$gene_family = factor(gen$gene_family, levels = unique(gen_list$Family))

# add species and taxonomic info, and factor order
gen = merge(gen, taxe, by.x = "species", by.y = "Species", all.x = T)
gen$Group_factor = factor(gen$Group, levels = unique(taxe$Group))
gen$Species_factor = factor(gen$Species.name, levels = unique(taxe$Species.name))
sps_gap_ixs = c(1,1+which(diff(as.numeric(taxe$Group))!=0)) - 1

# turn into a gene presence table (one entry per species)
gen_u = within(gen, rm("gene"))
gen_u = unique(gen_u)

# count number of species per taxonomic group
tax_group_counts = table(factor(taxe$Group, levels=unique(taxe$Group)))

# plot counts per group
pdf(file=sprintf("%s/counts_euk_genes-per_lineage.pdf", outp_fn),height=6,width=8)
gen_crosstab = xtabs(formula = ~ gene_family + Group_factor, data = gen_u, drop.unused.levels	= F)
pheatmap(t(gen_crosstab), color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         gaps_col = gen_list_gap_ixs,
         cellwidth = 5, cellheight = 5, na_col = "dodgerblue4", number_color = "aliceblue", fontsize = 5,
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_format = "%i",
         main=paste("Gene presence: Eukaryota lineages (num species)"))


# turn into fraction
gen_crosstab_frac = sweep(gen_crosstab, MARGIN=2, tax_group_counts, "/")
pheatmap(t(gen_crosstab_frac), 
         color = col_blue(20), breaks = c(0,0.0001,seq(0.05,1,length.out = 19)), 
         gaps_col = gen_list_gap_ixs,
         cellwidth = 5, cellheight = 5, na_col = "dodgerblue4", number_color = "aliceblue", fontsize = 5,
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = F, number_format = "%.2f",
         main=paste("Gene presence: Eukaryota lineages (fraction of species)"))

gen_crosstab = xtabs(formula = ~ gene_family + Group_factor, data = gen, drop.unused.levels	= F)
pheatmap(t(gen_crosstab), color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         gaps_col = gen_list_gap_ixs,
         cellwidth = 5, cellheight = 5, na_col = "dodgerblue4", number_color = "aliceblue", fontsize = 5,
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_format = "%i",
         main=paste("Gene presence: Eukaryota lineages (total num genes)"))

# plot num sps per group
par(mar=c(16, 4.1, 4.1, 2.1))
barplot(tax_group_counts, las=2, border = NA, names.arg = paste(names(tax_group_counts),"n =", tax_group_counts),
        main="num species per group")


# save csvs
write.table(t(gen_crosstab), file=sprintf("%s/counts_euk_genes_groups.csv",outp_fn), sep="\t", quote = F)
dev.off()

pdf(file=sprintf("%s/counts_euk_genes-per_species.pdf", outp_fn),height=16,width=8)

# plot counts per species
gen_crosstab = xtabs(formula = ~ gene_family + Species_factor, data = gen, drop.unused.levels	= F)
pheatmap(t(gen_crosstab), color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         gaps_col = gen_list_gap_ixs, gaps_row=sps_gap_ixs,
         cellwidth = 5, cellheight = 5, na_col = "dodgerblue4", number_color = "aliceblue", fontsize = 5,
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T, number_format = "%i",
         main=paste("Number of genes per species"))

# save csvs
write.table(t(gen_crosstab), file=sprintf("%s/counts_euk_genes_species.csv", outp_fn), sep="\t", quote = F)

dev.off()



#### eukaryotes OG counts per species ####

pdf(file=sprintf("%s/counts_euk_orthogroups-per_species.pdf", outp_fn),height=8,width=16)
gen_euktab = data.frame()
for (gei in as.vector(gen_list$Family)) {
    
    geit = data.frame()
    geit_list = Sys.glob(sprintf("%s/euk.%s.HG*.ortholog_groups.csv",phylo_fn, gei))
    if (length(geit_list) > 0) {
        
        print(sprintf("tables %s species",gei))
        
        # load all tables for this gene family
        for (geif in geit_list){
            geit = rbind(geit, read.table(geif, header = T))
        }
        geit = geit[  order( as.numeric(geit$orthogroup )),  ]
        geit$species = stringr::str_split(geit$gene, pattern = "_", simplify = T)[,1]
        geit$species_factor = factor(geit$species, levels = taxe$Species)
        geit = merge(geit, taxe, by.x = "species", by.y = "Species", all.x = T, all.y=F)
        geit$group_factor = factor(geit$Group, levels = unique(taxe$Group))
        
        # one row per species
        geit_u = within(geit, rm("gene"))
        geit_u = unique(geit_u)
        
        
        # cross tabulation of species and OGs
        gen_crosssps = xtabs(formula = ~ orthogroup + species_factor, data = geit, drop.unused.levels	= F)
        keep_non_singletons = rowSums(gen_crosssps>1) > 1
        gen_crosssps = matrix(gen_crosssps[keep_non_singletons,], nrow = length(which(keep_non_singletons)))
        rownames(gen_crosssps) = names(which(keep_non_singletons))
        colnames(gen_crosssps) = taxe$Species
        
        # plot
		gen_crosssps = order_matrix(gen_crosssps)
        pheatmap(gen_crosssps, 
                 color = col_blue(10), breaks = seq(0, 10, length.out = 11)-0.01, 
                 labels_col = stringr::str_trunc(taxe$Species.name, width = 50),
                 labels_row = stringr::str_trunc(rownames(gen_crosssps), width = 50),
                 gaps_col=sps_gap_ixs, legend = F,
                 cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5, 
                 border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T, number_format = "%i",
                 main=sprintf("%s OG presence", gei))
        
        # accumulate counts
        gen_euktab = rbind(gen_euktab, gen_crosssps)
        
    } else {
        
        print(sprintf("tables %s are empty, skip",gei))
    }
    
}

dev.off()
write.table(gen_euktab, file=sprintf("%s/counts_euk_orthogroups.csv", outp_fn), sep="\t", quote = F)



#### eukaryotes OG counts per phylum: OG presence counts ####

pdf(file=sprintf("%s/counts_euk_orthogroups-per_lineage-counts.pdf", outp_fn),height=6,width=8)
for (gei in as.vector(gen_list$Family)) {
    
    geit = data.frame()
    geit_list = Sys.glob(sprintf("%s/euk.%s.HG*.ortholog_groups.csv",phylo_fn, gei))
    if (length(geit_list) > 0) {
        
        print(sprintf("tables %s phyla",gei))
        
        # load all tables for this gene family
        for (geif in geit_list){
            geit = rbind(geit, read.table(geif, header = T))
        }
        geit = geit[  order( as.numeric(geit$orthogroup )),  ]
        geit$species = stringr::str_split(geit$gene, pattern = "_", simplify = T)[,1]
        geit$species_factor = factor(geit$species, levels = taxe$Species)
        geit = merge(geit, taxe, by.x = "species", by.y = "Species", all.x = T, all.y=F)
        geit$group_factor = factor(geit$Group, levels = unique(taxe$Group))
        
        # one row per species
        geit_u = within(geit, rm("gene"))
        geit_u = unique(geit_u)
        
        # find non-singletons at sps level
        gen_crosssps = xtabs(formula = ~ orthogroup + species_factor, data = geit, drop.unused.levels	= F)
        keep_non_singletons = rowSums(gen_crosssps>1) > 1
        
        # cross tabulation of taxon groups and OGs
        gen_crossgru = xtabs(formula = ~ orthogroup + group_factor, data = geit_u, drop.unused.levels	= F)
        gen_crossgru = matrix(gen_crossgru[keep_non_singletons,], nrow = length(which(keep_non_singletons)))
        rownames(gen_crossgru) = names(which(keep_non_singletons))
        colnames(gen_crossgru) = unique(taxe$Group)
        
        # plot
		gen_crossgru = order_matrix(gen_crossgru)
        pheatmap(t(gen_crossgru), color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
                 labels_row = unique(taxe$Group), legend = F,
                 labels_col = stringr::str_trunc(rownames(gen_crossgru), width = 50),
                 cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5, 
                 border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T, number_format = "%i",
                 main=sprintf("%s OG presence", gei))
        
        
    } else {
        
        print(sprintf("tables %s are empty, skip",gei))
    }
    
}
dev.off()


#### eukaryotes OG counts per phylum: OG presence fraction ####

pdf(file=sprintf("%s/counts_euk_orthogroups-per_lineage-frac.pdf", outp_fn),height=6,width=8)
for (gei in as.vector(gen_list$Family)) {
  
  geit = data.frame()
  geit_list = Sys.glob(sprintf("%s/euk.%s.HG*.ortholog_groups.csv",phylo_fn, gei))
  if (length(geit_list) > 0) {
    
    print(sprintf("tables %s phyla",gei))
    
    # load all tables for this gene family
    for (geif in geit_list){
      geit = rbind(geit, read.table(geif, header = T))
    }
    geit = geit[  order( as.numeric(geit$orthogroup )),  ]
    geit$species = stringr::str_split(geit$gene, pattern = "_", simplify = T)[,1]
    geit$species_factor = factor(geit$species, levels = taxe$Species)
    geit = merge(geit, taxe, by.x = "species", by.y = "Species", all.x = T, all.y=F)
    geit$group_factor = factor(geit$Group, levels = unique(taxe$Group))
    
    # one row per species
    geit_u = within(geit, rm("gene"))
    geit_u = unique(geit_u)
    
    # find non-singletons at sps level
    gen_crosssps = xtabs(formula = ~ orthogroup + species_factor, data = geit, drop.unused.levels	= F)
    keep_non_singletons = rowSums(gen_crosssps>1) > 1
    
    # cross tabulation of taxon groups and OGs
    gen_crossgru = xtabs(formula = ~ orthogroup + group_factor, data = geit_u, drop.unused.levels	= F)
    gen_crossgru_frac = sweep(gen_crossgru, MARGIN=1, tax_group_counts, "/")
    gen_crossgru_frac = matrix(gen_crossgru_frac[keep_non_singletons,], nrow = length(which(keep_non_singletons)))
    rownames(gen_crossgru_frac) = names(which(keep_non_singletons))
    colnames(gen_crossgru_frac) = unique(taxe$Group)
    
    # plot
    gen_crossgru_frac = order_matrix(gen_crossgru_frac)
    pheatmap(t(gen_crossgru_frac), 
             color = col_blue(20), breaks = c(0,0.0001,seq(0.05,1,length.out = 19)), 
             labels_row = unique(taxe$Group),legend = F,
             labels_col = stringr::str_trunc(rownames(gen_crossgru_frac), width = 50),
             cellwidth = 5, cellheight = 5, na_col = "dodgerblue4",number_color = "aliceblue", fontsize = 5, 
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F, number_format = "%i",
             main=sprintf("%s OG presence (fraction of species)", gei))
    
    
  } else {
    
    print(sprintf("tables %s are empty, skip",gei))
  }
  
}
dev.off()
