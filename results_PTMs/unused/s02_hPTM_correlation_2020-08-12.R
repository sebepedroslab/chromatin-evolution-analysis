# functions and libraries
library(pheatmap)
library(readxl)
library(gtools)
library(stringr)

graphics.off()

# input files
mod_fn = "consensus_modifications_perseq.xlsx"
sps_fn = "species_list_2020-07-22.txt"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
his_list = c("H2A","H2B","H3","H4","macroH2A","H2AZ")
gen_fn = "../results_toolkit_phylo/counts_orthogroups_eukaryotes.csv"

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
library(glmnet)
source("../scripts/summarise_model_OR.R")
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
  
  # empty matrix for correlation tests (pval)
  mod_cop = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  colnames(mod_cop) = mod_list
  rownames(mod_cop) = rownames(geb_m)
  
  # correlation tests for each modification
  for (moi in mod_list) {
    mod_cor[,moi] = apply(geb_m, 1, function(x) cor.test( x, unlist(mod[moi,]), method = "pearson")$estimate )
    mod_cop[,moi] = apply(geb_m, 1, function(x) cor.test( x, unlist(mod[moi,]), method = "pearson")$p.value )
  }
  
  # remove non-significant values from cor matrix
  mod_cor_net = mod_cor
  mod_cor_net [ mod_cop > 0.05 ] = NA
  
  mod_cor_rowclus = hclust(dist(mod_cor), method = "ward.D2", )
  # plot correlation tests
  pdf(sprintf("models_Pearson.Ac.%s.pdf", his),height=6,width=8)
  pheatmap(mod_cor, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, Pearson correlation", his))
  
  # plot correlation tests (no non-sig)
  pheatmap(mod_cor_net, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, Pearson correlation", his))
  
  # plot some sort of overlap matrix?
  for (moi in mod_list) {
    
    cor_order = order(mod_cor[,moi], decreasing = T)
    mod_ove = rbind(mod[moi,],geb_m[cor_order,])
    annotation_df = rbind(-1,data.frame(mod_cor_net[,moi]))
    rownames(annotation_df)[1] = moi
    pheatmap(mod_ove, color = cbw.fun(2), breaks = seq(0,1,length.out = 3), 
             gaps_row = 1, annotation_row = annotation_df,
             cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
             main=sprintf("%s %s hPTMs and enzymes", his, moi))
    
  }
  dev.off()
  
}


#### 3.2 Correlation: methylation ####

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
  
  # empty matrix for correlation tests (pval)
  mod_cop = matrix(nrow=nrow(geb_m), ncol = length(mod_list))
  colnames(mod_cop) = mod_list
  rownames(mod_cop) = rownames(geb_m)
  
  # correlation tests for each modification
  for (moi in mod_list) {
    mod_cor[,moi] = apply(geb_m, 1, function(x) cor.test( x, unlist(mod[moi,]), method = "pearson")$estimate )
    mod_cop[,moi] = apply(geb_m, 1, function(x) cor.test( x, unlist(mod[moi,]), method = "pearson")$p.value )
  }
  
  # remove non-significant values from cor matrix
  mod_cor_net = mod_cor
  mod_cor_net [ mod_cop > 0.05 ] = NA
  
  mod_cor_rowclus = hclust(dist(mod_cor), method = "ward.D2", )
  # plot correlation tests
  pdf(sprintf("models_Pearson.Me.%s.pdf", his),height=6,width=8)
  pheatmap(mod_cor, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, Pearson correlation", his))
  
  # plot correlation tests (no non-sig)
  pheatmap(mod_cor_net, color = cod.fun(31), breaks = seq(-1,1,length.out = 32), 
           cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
           border_color = "white", cluster_cols=F, cluster_rows=mod_cor_rowclus, display_numbers = F,
           main=sprintf("%s hPTMs, Pearson correlation", his))
  
  # plot some sort of overlap matrix?
  for (moi in mod_list) {
    
    cor_order = order(mod_cor[,moi], decreasing = T)
    mod_ove = rbind(mod[moi,],geb_m[cor_order,])
    annotation_df = rbind(-1,data.frame(mod_cor_net[,moi]))
    rownames(annotation_df)[1] = moi
    pheatmap(mod_ove, color = cbw.fun(2), breaks = seq(0,1,length.out = 3), 
             gaps_row = 1, annotation_row = annotation_df,
             cellwidth = 6, cellheight = 6, na_col = "grey", fontsize = 5,
             border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
             main=sprintf("%s %s hPTMs and enzymes", his, moi))
    
  }
  dev.off()
  
}

#### 3.1 Acetylation regression ####

pdf(sprintf("models_LASSO.Ac.pdf", his),height=6,width=8)

for (his in his_list) {
  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  mod = mod[grep("Ace$", rownames(mod)),]
  mod_list = rownames(mod)
  print(sprintf("Report %s correlation", his))
  
  # loop through marks and apply LASSO variant selection
  # output matrix
  coef_matrix_lasso = matrix(nrow = nrow(geb_ac), ncol = length(mod_list))
  rownames(coef_matrix_lasso) = rownames(geb_ac)
  colnames(coef_matrix_lasso) = mod_list
  
  for (moi in mod_list) {
    
    # data
    mark = mod[moi,]
    is_mark_na = is.na(mark)
    mark = mark[!is_mark_na]
    data = t(geb_ac[,!is_mark_na])
    
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
  # coef_matrix_lasso_df = coef_matrix_lasso_df[rowSums(coef_matrix_lasso_df)>0,]
  # coef_matrix_lasso_df = coef_matrix_lasso
  
  # print coeficients in a heatmap
  if (ncol(coef_matrix_lasso_df) > 0) {
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


#### 3.2 Methylation regression ####

pdf(sprintf("models_LASSO.Me.pdf", his),height=6,width=8)

for (his in his_list) {
  
  mod = read.table(sprintf("modifications_per_sps_%s.csv", his), header = T, row.names = 1, sep = "\t")
  mod = mod[grep("Me[1-3]$", rownames(mod)),]
  mod_list = rownames(mod)
  
  # loop through marks and apply LASSO variant selection
  # output matrix
  coef_matrix_lasso = matrix(nrow = nrow(geb_me), ncol = length(mod_list))
  rownames(coef_matrix_lasso) = rownames(geb_me)
  colnames(coef_matrix_lasso) = mod_list
  
  for (moi in mod_list) {
    
    # data
    mark = mod[moi,]
    is_mark_na = is.na(mark)
    mark = mark[!is_mark_na]
    data = t(geb_me[,!is_mark_na])
    
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
  # coef_matrix_lasso_df = coef_matrix_lasso_df[rowSums(coef_matrix_lasso_df)>0,]
  # coef_matrix_lasso_df = coef_matrix_lasso
  
  # print coeficients in a heatmap
  if (ncol(coef_matrix_lasso_df) > 0) {
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


stop("ARA")


#### 2.1 Acetylation regression ####

mod = mod_ac
mod_list = rownames(mod_ac)

# loop through marks and apply BIC variant removal
# output matrix
coef_matrix_BIC_ac = matrix(nrow = nrow(geb_ac), ncol = length(mod_list))
rownames(coef_matrix_BIC_ac) = rownames(geb_ac)
colnames(coef_matrix_BIC_ac) = mod_list

for (moi in mod_list) {
  
  data = as.data.frame(t(geb_ac))
  data$mark = unlist(mod[moi,])
  data = as.data.frame(data)
  
  if (length(unique(data[,"mark"])) > 1) {
    
    print(sprintf("%s regression analysis", moi))
    
    # full model
    mod_tot = glm(mark ~ ., data = data, family = "binomial")
    # reduce full model using BIC
    mod_min = step(mod_tot, direction = "both", steps = 1000, trace = F, k = log(nrow(data)), ) # k=log(num_obs) for BIC
    # null model
    mod_nul = glm(mark ~ 1, data = data, family = "binomial")
    
    # report min BIC model
    mod_min_signif = glm_tables(
      model = mod_min,
      null = mod_nul, 
      model_name = moi)
    
    write.table(file="models_BIC.Ac.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="models_BIC.Ac.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    write.table(file="models_BIC.Ac.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
    # model coefficients
    model_min_coefs = as.matrix(coef(mod_min))
    model_min_coefs_sig_names = rownames(model_min_coefs)[-1]
    
    # add coefs to final matrix
    coef_matrix_BIC_ac[gsub("`","",model_min_coefs_sig_names),moi] = model_min_coefs[-1]
  } else {
    print(sprintf("%s is constant, skip regression", moi))
    print(mod[moi,])
  }
  
}

# print coeficients in a heatmap
pdf("models_BIC.Ac.pdf",height=8,width=16)
coef_matrix_lasso_df = coef_matrix_BIC_ac[rowSums(is.na(coef_matrix_BIC_ac))<ncol(coef_matrix_BIC_ac),colSums(is.na(coef_matrix_BIC_ac))<nrow(coef_matrix_BIC_ac)]
pheatmap(coef_matrix_lasso_df, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
         cellwidth = 6, cellheight = 6, na_col = "grey",number_color = "green", fontsize = 6,
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
         main="coefficients BIC regression, Ace")
dev.off()


#### 2.2 Methylation regression ####

mod_list = rownames(mod_me)

# loop through marks and apply BIC variant removal
# output matrix
coef_matrix_BIC_me = matrix(nrow = nrow(geb_me), ncol = length(mod_list))
rownames(coef_matrix_BIC_me) = rownames(geb_me)
colnames(coef_matrix_BIC_me) = mod_list

for (moi in mod_list) {
  
  data = as.data.frame(t(geb_me))
  data$mark = unlist(mod_me[moi,])
  data = as.data.frame(data)
  
  if (length(unique(data[,"mark"])) > 1) {
    
    print(sprintf("%s regression analysis", moi))
    
    # full model
    mod_tot = glm(mark ~ ., data = data, family = "binomial")
    # reduce full model using BIC
    mod_min = step(mod_tot, direction = "both", steps = 1000, trace = F, k = log(nrow(data)), ) # k=log(num_obs) for BIC
    # null model
    mod_nul = glm(mark ~ 1, data = data, family = "binomial")
    
    # report min BIC model
    mod_min_signif = glm_tables(
      model = mod_min,
      null = mod_nul, 
      model_name = moi)
    
    write.table(file="models_BIC.Me.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="models_BIC.Me.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    write.table(file="models_BIC.Me.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
    # model coefficients
    model_min_coefs = as.matrix(coef(mod_min))
    model_min_coefs_sig_names = rownames(model_min_coefs)[-1]
    
    # add coefs to final matrix
    coef_matrix_BIC_me[gsub("`","",model_min_coefs_sig_names),moi] = model_min_coefs[-1]
  } else {
    print(sprintf("%s is constant, skip regression", moi))
    print(mod_me[moi,])
  }
  
}

# print coeficients in a heatmap
pdf("~/models_BIC.Me.pdf",height=8,width=16)
coef_matrix_lasso_df = coef_matrix_BIC_me[rowSums(is.na(coef_matrix_BIC_me))<ncol(coef_matrix_BIC_me),colSums(is.na(coef_matrix_BIC_me))<nrow(coef_matrix_BIC_me)]
pheatmap(coef_matrix_lasso_df, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
         cellwidth = 6, cellheight = 6, na_col = "grey",number_color = "green", fontsize = 6,
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F,
         main="coefficients BIC regression, Met")
dev.off()





stop("AERA")

