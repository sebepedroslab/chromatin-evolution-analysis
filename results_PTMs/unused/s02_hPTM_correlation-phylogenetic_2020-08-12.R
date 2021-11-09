# functions and libraries
library(pheatmap)
library(readxl)
library(gtools)
library(stringr)
library(phyr)
library(ape)

graphics.off()

# input files
mod_fn = "consensus_modifications_perseq.xlsx"
sps_fn = "species_list_2020-07-22.txt"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
his_list = c("H2A","H2B","H3","H4","macroH2A","H2AZ")
gen_fn = "../results_toolkit_phylo/counts_orthogroups_eukaryotes.csv"
phyl_fn = "species_list_2020-07-22.newick"

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


#### 2. Regression analysis ####


mod_me = data.frame()
mod_ac = data.frame()
for (his in his_list) {
  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  
  # get methylation events
  mod_mi = mod[grep("Me[1-3]$", rownames(mod)),]
  rownames(mod_mi) = paste(his, rownames(mod_mi))
  mod_me = rbind(mod_me, mod_mi)
  
  # get acetylation events
  mod_ai = mod[grep("Ace$", rownames(mod)),]
  rownames(mod_ai) = paste(his, rownames(mod_ai))
  mod_ac = rbind(mod_ac, mod_ai)
  
  
}

# functions and libraries
library(pheatmap)
# library(glmnet)
# source("../scripts/summarise_model_OR.R")
library(plyr)


# input files: all orthogroups
gen = read.table(gen_fn, header = T, sep = "\t")
# remove species without hPTM data
geb = gen[,colnames(gen) %in% sps_list_ab]
# binarise orthogroup presence
geb = (geb>1) * 1
# remove orthogroups that are always absent, and remove singletons
geb = geb[rowSums(geb)>0,]
geb = geb[rowSums(geb)>1,]

# load list hmms
gen_list_fn = "../data/gene_families_hmm.csv"
gen_list = read.table(gen_list_fn, header = F)
colnames(gen_list) = c("Class","Type","Family","Domains","search_thr","inflation","min_phylo_size")

# get separate tables for acetylation enzymes
gen_list_ac_doms = as.vector(gen_list[gen_list$Type == "Acetylase","Domains"])
gen_list_ac_doms = as.vector(gen_list[gen_list$Class == "Acetylation","Domains"])
gen_list_ac_doms = unlist(stringr::str_split(string = gen_list_ac_doms, pattern = ","))
geb_ac = geb[grep(pattern = paste(gen_list_ac_doms, collapse = "|"), x = rownames(geb)),]

# get separate tables for methylation enzymes
gen_list_me_doms = as.vector(gen_list[gen_list$Type == "Methylase","Domains"])
gen_list_me_doms = as.vector(gen_list[gen_list$Class == "Methylation","Domains"])
gen_list_me_doms = unlist(stringr::str_split(string = gen_list_me_doms, pattern = ","))
geb_me = geb[grep(pattern = paste(gen_list_me_doms, collapse = "|"), x = rownames(geb)),]


# read phylogeny
phyl = ape::read.tree(phyl_fn)


#### 2. Correlation ####

#### 2.1 Correlation: acetylation ####

for (his in his_list) {
  
  # load data for this particular histone  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  mod = mod[grep("Ace$", rownames(mod)),]
  mod_list = rownames(mod)
  geb_m = geb_ac
  print(sprintf("Report %s correlation", his))
  
  
  # empty matrix for correlation tests
  mod_cor = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  colnames(mod_cor) = mod_list
  rownames(mod_cor) = rownames(geb_m)
  
  # # empty matrix for correlation tests (pval)
  # mod_cop = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  # colnames(mod_cop) = mod_list
  # rownames(mod_cop) = rownames(geb_m)

  for (moi in mod_list) {
    for (gei in rownames(geb_m)) {
      geb_mi = data.frame(geb_m[gei,])
      geb_mi$mod = unlist(mod[moi,])
      colnames(geb_mi) = c("gene","mod")
      formula_string = paste("~",paste(colnames(geb_mi), collapse = " + "))
      if (sd(geb_mi$mod>0)) { 
        geb_mi_corphy = cor_phylo(
          variates = as.formula(formula_string),
          species = rownames(geb_mi), 
          data = geb_mi,
          phy = phyl)
        # store result
        mod_cor[gei,moi] = geb_mi_corphy$corrs[1,2]
        # mod_cop[gei,moi] = geb_mi_corphy$B[1,"P-value"]
      } else { 
        mod_cor[gei,moi] = NA
        # mod_cop[gei,moi] = NA
      }
    }
  }
  
  # # remove non-significant values from cor matrix
  # mod_cor_net = mod_cor
  # mod_cor_net [ mod_cop >= 0.05 ] = NA
  
  mod_cor_rowclus = hclust(dist(mod_cor), method = "ward.D2")
  # plot correlation tests
  pdf(sprintf("models_phyloPearson.Ac.%s.pdf", his),height=6,width=8)
  pheatmap(mod_cor, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, phylogenetic Pearson correlation", his))
  
  # pheatmap(mod_cor_net, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
  #          cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
  #          border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
  #          main=sprintf("%s hPTMs, phylogenetic Pearson correlation", his))
  
  # plot some sort of overlap matrix?
  for (moi in mod_list) {
    
    cor_order = order(mod_cor[,moi], decreasing = T)
    mod_ove = rbind(mod[moi,],geb_m[cor_order,])
    annotation_df = rbind(-1,data.frame(mod_cor[,moi]))
    rownames(annotation_df)[1] = moi
    pheatmap(mod_ove, color = cbw.fun(2), breaks = seq(0,1,length.out = 3), 
             gaps_row = 1, annotation_row = annotation_df,
             cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
             main=sprintf("%s %s hPTMs and enzymes", his, moi))
    
  }
  dev.off()
  
}


#### 2.2 Correlation: methylation ####

for (his in his_list) {
  
  # load data for this particular histone  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  mod = mod[grep("Me[1-3]$", rownames(mod)),]
  mod_list = rownames(mod)
  geb_m = geb_me
  print(sprintf("Report %s correlation", his))
  
  
  # empty matrix for correlation tests
  mod_cor = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  colnames(mod_cor) = mod_list
  rownames(mod_cor) = rownames(geb_m)
  
  # # empty matrix for correlation tests (pval)
  # mod_cop = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  # colnames(mod_cop) = mod_list
  # rownames(mod_cop) = rownames(geb_m)
  
  for (moi in mod_list) {
    for (gei in rownames(geb_m)) {
      geb_mi = data.frame(geb_m[gei,])
      geb_mi$mod = unlist(mod[moi,])
      colnames(geb_mi) = c("gene","mod")
      formula_string = paste("~",paste(colnames(geb_mi), collapse = " + "))
      if (sd(geb_mi$mod>0)) { 
        geb_mi_corphy = cor_phylo(
          variates = as.formula(formula_string),
          species = rownames(geb_mi), 
          data = geb_mi,
          phy = phyl)
        # store result
        mod_cor[gei,moi] = geb_mi_corphy$corrs[1,2]
        # mod_cop[gei,moi] = geb_mi_corphy$B[1,"P-value"]
      } else { 
        mod_cor[gei,moi] = NA
        # mod_cop[gei,moi] = NA
      }
    }
  }
  
  # # remove non-significant values from cor matrix
  # mod_cor_net = mod_cor
  # mod_cor_net [ mod_cop >= 0.05 ] = NA
  
  mod_cor_rowclus = hclust(dist(mod_cor), method = "ward.D2")
  # plot correlation tests
  pdf(sprintf("models_phyloPearson.Me.%s.pdf", his),height=6,width=8)
  pheatmap(mod_cor, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, phylogenetic Pearson correlation", his))
  
  # pheatmap(mod_cor_net, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
  #          cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
  #          border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
  #          main=sprintf("%s hPTMs, phylogenetic Pearson correlation", his))
  
  # plot some sort of overlap matrix?
  for (moi in mod_list) {
    
    cor_order = order(mod_cor[,moi], decreasing = T)
    mod_ove = rbind(mod[moi,],geb_m[cor_order,])
    annotation_df = rbind(-1,data.frame(mod_cor[,moi]))
    rownames(annotation_df)[1] = moi
    pheatmap(mod_ove, color = cbw.fun(2), breaks = seq(0,1,length.out = 3), 
             gaps_row = 1, annotation_row = annotation_df,
             cellwidth = 5, cellheight = 5, na_col = "grey", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
             main=sprintf("%s %s hPTMs and enzymes", his, moi))
    
  }
  dev.off()
  
}


