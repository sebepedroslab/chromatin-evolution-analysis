# functions and libraries
library(pheatmap)
library(readxl)
library(gtools)
library(stringr)

graphics.off()

# input files: Methylation
mod_fn = "consensus_modifications_perseq.xlsx"
sps_fn = "species_list_2020-07-22.txt"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
his_list = c("H2A","H2B","H3","H4","macroH2A","H2AZ")

# read ref seqs
mod_ref = read_excel(mod_fn, sheet="ref_seqs")

# heatmap color
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
col.fun = colorRampPalette(interpolate="l",c("azure2","deepskyblue","dodgerblue4"))
cod.fun = colorRampPalette(interpolate="l",c("firebrick4","firebrick3","orangered", "gray95", "deepskyblue","dodgerblue3","dodgerblue4"))


# read sps list (important for sps order)
sps = read.table(sps_fn, stringsAsFactors = F, header = F, col.names = c("species"), sep = "\t")
sps_order = sps$species
tax = read.table(tax_fn, stringsAsFactors = F, sep = "\t", header = T)
sps = merge(sps, tax, by.x="species", by.y="Species", all.x=T)
sps = sps[match(sps_order, sps$species),]
sps_list = sps$Species.name
sps_list_ab = sps$species




#### 1. Presence/absence of hPTMs ####

pdf("modifications_per_sps.pdf",height=6,width=16)
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
  
  # reorder and get modification names
  mod = mod[  order( as.numeric(mod$Positions_in_Reference) , mod$Positions_in_Reference , mod$Modification ),  ]
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
  
  
  keep_nonsingletons = rowSums(mod_per_sps>1)>1
  
  
  # plot
  pheatmap(t(mod_per_sps), color = col.fun(2), breaks = seq(0,1,length.out = 3), 
           gaps_col = mod_per_sps_pos_gaps, gaps_row = 13, fontsize = 5,
           cellwidth = 5, cellheight = 5, na_col = "grey", number_color = "aliceblue", number_format = "%i",
           border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T,
           main=sprintf("Histone %s", his))
  
  # accumulate counts
  mod_per_sps_m = matrix(mod_per_sps, nrow = nrow(mod_per_sps))
  mod_per_sps_m = (mod_per_sps_m > 0) * 1 # binarise
  colnames(mod_per_sps_m) = sps_list_ab
  rownames(mod_per_sps_m) = stringr::str_split(rownames(mod_per_sps), pattern = " \\(", simplify = T)[,1]
  write.table(mod_per_sps_m, file=sprintf("modifications_per_sps_%s.csv", his), sep="\t", quote = F)
  
}
dev.off()



#### 2. Regression analysis ####

# functions and libraries
library(pheatmap)
library(glmnet)
source("../scripts/summarise_model_OR.R")
library(plyr)

# input files: all orthogroups
gen_fn = "../results_toolkit_phylo/counts_orthogroups_eukaryotes.csv"
gen = read.table(gen_fn, header = T, sep = "\t")

# remove species without hPTM data
geb = gen[,colnames(gen) %in% sps_list_ab]

pdf(sprintf("models_LASSO.pdf", his),height=6,width=8)
for (his in his_list) {
  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  mod_list = rownames(mod)
  
  # loop through marks and apply LASSO variant selection
  # output matrix
  coef_matrix_lasso = matrix(nrow = nrow(geb), ncol = length(mod_list))
  rownames(coef_matrix_lasso) = rownames(geb)
  colnames(coef_matrix_lasso) = mod_list

  for (moi in mod_list) {
    
    # data
    mark = mod[moi,]
    is_mark_na = is.na(mark)
    mark = mark[!is_mark_na]
    data = t(geb[,!is_mark_na])
    
    if (length(unique(mark)) > 1 && min(table(mark))>3) {
      
      # cross-validated model to find best lambda for lasso
      model_cv = cv.glmnet(x=data, y=mark, family='binomial')
      
      # recalculate model as glmnet object, with min lambda
      model_lasso = glmnet(x=data, y=mark, family = "binomial", lambda = model_cv$lambda.min, alpha = 1)
      
      # find significant variants (non-zero coefficients)
      model_lasso_coefs = as.matrix(coef(model_lasso))
      model_lasso_coefs_sig = model_lasso_coefs[model_lasso_coefs!=0,]
      model_lasso_coefs_sig_names = names(model_lasso_coefs_sig)[-1]
      data_min = as.data.frame(data[,model_lasso_coefs_sig_names])
      data_min[,"mark"] = mark
      
      # dummy min model with only minimal data
      model_lasso_dummy = glm(mark ~., data = data_min, family = "binomial")
      
      # null model
      model_null = glm(mark ~ 1, data = data_min, family = "binomial")
      
      # model significance
      mod_min_signif = glm_tables(
        model=model_lasso_dummy, 
        null = model_null,
        model_name = moi)
      
      # # report model
      # write.table(file="models_LASSO.Me.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
      # write.table(file="models_LASSO.Me.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
      # write.table(file="models_LASSO.Me.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
      
      # add coefs to final matrix
      coef_matrix_lasso[,moi] = model_lasso_coefs[-1]
    }
  }
  
  coef_matrix_lasso_df = coef_matrix_lasso[,colSums(is.na(coef_matrix_lasso))<nrow(coef_matrix_lasso)]
  coef_matrix_lasso_df = coef_matrix_lasso_df[rowSums(coef_matrix_lasso_df)>0,]
    
  # print coeficients in a heatmap
  if (nrow(coef_matrix_lasso_df > 0)) {
    print(sprintf("Report %s regression analysis", his))
    pheatmap(coef_matrix_lasso_df, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
           cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
           main=sprintf("%s hPTMs, LASSO coefficients", his))
  } else {
    print(sprintf("Skip %s regression analysis", his))
  }
    
}
dev.off()


