# functions and libraries
library(pheatmap)
library(readxl)
library(gtools)
library(stringr)

graphics.off()

# input files
mod_fn = "consensus_modifications_perseq.xlsx"
sps_fn = "species_list_2021-02-10.txt"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
his_list = c("ubH3","ubH4","ubH2A","ubH2B","ubmacroH2A","ubH2AZ","ubH1")

# read ref seqs
mod_ref = read_excel(mod_fn, sheet="ref_seqs")

# heatmap color
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
col.fun = colorRampPalette(interpolate="l",c("azure2","deepskyblue","dodgerblue4"))
cod.fun = colorRampPalette(interpolate="l",c("firebrick4","firebrick3","orangered", "gray95", "deepskyblue","dodgerblue3","dodgerblue4"))
cbw.fun = colorRampPalette(interpolate="l",c("gray90","gray10"))


# read sps list (important for sps order)
sps = read.table(sps_fn, stringsAsFactors = F, header = F, col.names = c("species"), sep = "\t")
sps_order = sps$species
tax = read.table(tax_fn, stringsAsFactors = F, sep = "\t", header = T)
sps = merge(sps, tax, by.x="species", by.y="Species", all.x=T)
sps = sps[match(sps_order, sps$species),]
sps_list = sps$Species.name
sps_list_ab = sps$species




#### 1. Presence/absence of hPTMs ####

pdf("ubiquitination_per_sps.pdf",height=4,width=5)
mod_pres = data.frame()
for (his in his_list) {
    
    print(sprintf("Report %s counts", his))
    
    refseq = as.character(mod_ref[mod_ref$histone==his,"ref_seq"])
    mod = read_excel(mod_fn, sheet=his)
    mod = mod[!is.na(mod$Positions_in_Reference),]
    mod = mod[is.na(mod$ptmRS_Best_Site_Probabilities) | mod$ptmRS_Best_Site_Probabilities > 45,]
    # mod$Positions_in_Reference = as.numeric(mod$Positions_in_Reference)
    
    # add species
    mod$species = stringr::str_split(mod$Protein_ID, pattern = "_", simplify = T)[,1]
    mod = merge(mod,sps, by.x = "species", by.y="species", all.x=T, all.y=F)
    mod$Modification = factor(mod$Modification, levels = c("Ace","Me1","Me2","Me3","Cro","Pho","Ubi"))
    
    # reorder and get modification names
    mod = mod[  order( as.numeric(mod$Positions_in_Reference), mod$Positions_in_Reference , mod$Modification ),  ]
    # mod$modification_id = paste(
    #   paste(mod$Residue, mod$Positions_in_Reference, sep=""), 
    #   mod$Modification,
    #   paste("(ref:",mod$aa_in_reference,")", sep=""))
    
    mod$modification_id = paste(
        paste(mod$Residue, mod$Positions_in_Reference, sep=""), 
        mod$Modification, sep=" ")
    
    # cross-tabulate mods and species
    mod$species_factor = factor(mod$Species.name,levels = sps_list)
    mod_per_sps = xtabs(formula =  ~ modification_id + species_factor, data = mod)
    
    # order mods
    order_mods = unique(mod$modification_id)
    mod_per_sps = mod_per_sps[match(order_mods, rownames(mod_per_sps)),]
    
    # add gaps to rows
    mod_per_sps_pos = gsub( "^[A-Za-z]+", "", str_split(rownames(mod_per_sps), pattern = " ", simplify = T)[,1] )
    mod_per_sps_pos = as.numeric(mod_per_sps_pos)
    mod_per_sps_pos_gaps = c(1,1+which(diff(mod_per_sps_pos)!=0)) - 1
    mod_per_sps_pos_gaps = c(mod_per_sps_pos_gaps, which(is.na(mod_per_sps_pos)) -1 )
    
    # plot
    pheatmap(t(mod_per_sps), color = col.fun(2), breaks = seq(0,1,length.out = 3), 
             gaps_col = mod_per_sps_pos_gaps, gaps_row = c(12,24), fontsize = 5,
             cellwidth = 5, cellheight = 5, na_col = "grey", number_color = "aliceblue", number_format = "%i",
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T,
             main=sprintf("Histone %s", his))
    
    # accumulate counts
    mod_per_sps_m = matrix(mod_per_sps, nrow = nrow(mod_per_sps))
    mod_per_sps_m = (mod_per_sps_m > 0) * 1 # binarise
    colnames(mod_per_sps_m) = sps_list_ab
    rownames(mod_per_sps_m) = stringr::str_split(rownames(mod_per_sps), pattern = " \\(", simplify = T)[,1]
    write.table(mod_per_sps_m, file=sprintf("ubiquitination_per_sps_%s.csv", his), sep="\t", quote = F)
    
}
dev.off()


pdf("ubiquitination_per_sps-conserved_only.pdf",height=4,width=4)
mod_pres = data.frame()
for (his in his_list) {
    
    print(sprintf("Report %s counts", his))
    
    refseq = as.character(mod_ref[mod_ref$histone==his,"ref_seq"])
    mod = read_excel(mod_fn, sheet=his)
    mod = mod[!is.na(mod$Positions_in_Reference),]
    # mod$Positions_in_Reference = as.numeric(mod$Positions_in_Reference)
    
    # add species
    mod$species = stringr::str_split(mod$Protein_ID, pattern = "_", simplify = T)[,1]
    mod = merge(mod,sps, by.x = "species", by.y="species", all.x=T, all.y=F)
    mod$Modification = factor(mod$Modification, levels = c("Ace","Me1","Me2","Me3","Cro","Pho","Ubi"))
    
    # reorder and get modification names
    mod = mod[  order( as.numeric(mod$Positions_in_Reference), mod$Positions_in_Reference , mod$Modification ),  ]
    # mod$modification_id = paste(
    #   paste(mod$Residue, mod$Positions_in_Reference, sep=""), 
    #   mod$Modification,
    #   paste("(ref:",mod$aa_in_reference,")", sep=""))
    
    mod$modification_id = paste(
        paste(mod$Residue, mod$Positions_in_Reference, sep=""), 
        mod$Modification, sep=" ")
    
    # cross-tabulate mods and species
    mod$species_factor = factor(mod$Species.name,levels = sps_list)
    mod_per_sps = xtabs(formula =  ~ modification_id + species_factor, data = mod)
    
    # order mods
    order_mods = unique(mod$modification_id)
    mod_per_sps = mod_per_sps[match(order_mods, rownames(mod_per_sps)),]
    
    keep_nonsingletons = rowSums(mod_per_sps>0)>1
    keep_nonsingletons[is.na(keep_nonsingletons)] = F
    # list_residues_nonsing = paste(unique(stringr::str_split(names(keep_nonsingletons)[keep_nonsingletons], " ", simplify = T)[,1]), collapse = "|")
    # keep_nonsingletons_allmods = grepl(list_residues_nonsing, names(keep_nonsingletons))
    mod_per_sps = mod_per_sps[keep_nonsingletons,]
    
    if (sum(keep_nonsingletons)>1) {
        
        # add gaps to rows
        mod_per_sps_pos = gsub( "^[A-Za-z]+", "", str_split(rownames(mod_per_sps), pattern = " ", simplify = T)[,1] )
        mod_per_sps_pos = as.numeric(mod_per_sps_pos)
        mod_per_sps_pos_gaps = c(1,1+which(diff(mod_per_sps_pos)!=0)) - 1
        mod_per_sps_pos_gaps = c(mod_per_sps_pos_gaps, which(is.na(mod_per_sps_pos)) -1 )
        
        # plot
        pheatmap(t(mod_per_sps), color = col.fun(2), breaks = seq(0,1,length.out = 3), 
                 gaps_col = mod_per_sps_pos_gaps, gaps_row = c(12,24), fontsize = 4,
                 cellwidth = 4, cellheight = 4, na_col = "grey", number_color = "aliceblue", number_format = "%i",
                 border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
                 main=sprintf("Histone %s", his))
        
        # accumulate counts
        mod_per_sps_m = matrix(mod_per_sps, nrow = nrow(mod_per_sps))
        mod_per_sps_m = (mod_per_sps_m > 0) * 1 # binarise
        colnames(mod_per_sps_m) = sps_list_ab
        rownames(mod_per_sps_m) = stringr::str_split(rownames(mod_per_sps), pattern = " \\(", simplify = T)[,1]
        write.table(mod_per_sps_m, file=sprintf("ubiquitination_per_sps_%s.csv", his), sep="\t", quote = F)
    } else {
        print("skip, nothing is conserved...")
        
    }
    
}
dev.off()
