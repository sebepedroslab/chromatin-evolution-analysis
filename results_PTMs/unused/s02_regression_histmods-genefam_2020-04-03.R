# functions and libraries
library(pheatmap)
library(glmnet)
source("../scripts/summarise_model_OR.R")
library(plyr)

# input files: Methylation
mod_fn = ""
gen_fn = "../results_toolkit_phylo/counts_orthogroups_eukaryotes.csv"

# heatmap colors
col.fun = colorRampPalette(interpolate="l",c("azure","deepskyblue","dodgerblue4"))
cod.fun = colorRampPalette(interpolate="l",c("firebrick4","firebrick3","orangered", "gray95", "deepskyblue","dodgerblue3","dodgerblue4"))

# load data
gen = read.table(gen_fn, header = T, row.names = 1)
mod = read.table(mod_fn, header = T, row.names = 1)
sps_list = rownames(mod)
mod_list = colnames(mod)

# binarise gene presence
geb = gen>7
geb = geb * 1

# subset gene presence to species with hist mod info
geb = geb[sps_list,]

# loop through marks and apply LASSO variant selection
# output matrix
coef_matrix_lasso = matrix(nrow = ncol(geb), ncol = length(mod_list))
rownames(coef_matrix_lasso) = colnames(geb)
colnames(coef_matrix_lasso) = mod_list
for (moi in mod_list) {
  
  # data
  mark = mod[,moi]
  is_mark_na = is.na(mark)
  mark = mark[!is_mark_na]
  data = geb[!is_mark_na,]
  
  if (length(unique(mark)) > 1 && min(table(mark))>2) {
    
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
    
    # report model
    write.table(file="models_LASSO.Me.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="models_LASSO.Me.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    write.table(file="models_LASSO.Me.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
    # add coefs to final matrix
    coef_matrix_lasso[,moi] = model_lasso_coefs[-1]
  }
}

# print coeficients in a heatmap
pdf("models_LASSO.Me.pdf",height=20,width=20)
pheatmap(coef_matrix_lasso, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
         cellwidth = 16, cellheight = 16, na_col = "grey",number_color = "darkgreen",
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T,
         main="coefficients LASSO regression, Me")
dev.off()






# loop through marks and apply ridge variant selection
# output matrix
coef_matrix_ridge = matrix(nrow = ncol(geb), ncol = length(mod_list))
rownames(coef_matrix_ridge) = colnames(geb)
colnames(coef_matrix_ridge) = mod_list
for (moi in mod_list) {
  
  # data
  mark = mod[,moi]
  is_mark_na = is.na(mark)
  mark = mark[!is_mark_na]
  data = geb[!is_mark_na,]
  
  if (length(unique(mark)) > 1 && min(table(mark))>2) {
    
    # cross-validated model to find best lambda for ridge
    model_cv = cv.glmnet(x=data, y=mark, family='binomial')
    
    # recalculate model as glmnet object, with min lambda
    model_ridge = glmnet(x=data, y=mark, family = "binomial", lambda = model_cv$lambda.min, alpha = 0)
    
    # find significant variants (non-zero coefficients)
    model_ridge_coefs = as.matrix(coef(model_ridge))
    model_ridge_coefs_sig = model_ridge_coefs[model_ridge_coefs!=0,]
    model_ridge_coefs_sig_names = names(model_ridge_coefs_sig)[-1]
    data_min = as.data.frame(data[,model_ridge_coefs_sig_names])
    data_min[,"mark"] = mark
    
    # dummy min model with only minimal data
    model_ridge_dummy = glm(mark ~., data = data_min, family = "binomial")
    
    # null model
    model_null = glm(mark ~ 1, data = data_min, family = "binomial")
    
    # model significance
    mod_min_signif = glm_tables(
      model=model_ridge_dummy, 
      null = model_null,
      model_name = moi)
    
    # report model
    write.table(file="models_ridge.Me.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="models_ridge.Me.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    write.table(file="models_ridge.Me.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
    # add coefs to final matrix
    coef_matrix_ridge[,moi] = model_ridge_coefs[-1]
  }
}

# print coeficients in a heatmap
pdf("models_ridge.Me.pdf",height=20,width=20)
pheatmap(coef_matrix_ridge, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
         cellwidth = 16, cellheight = 16, na_col = "grey",number_color = "darkgreen",
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T,
         main="coefficients ridge regression, Me")
dev.off()




# loop through marks and apply BIC variant removal
# output matrix
coef_matrix_BIC = matrix(nrow = ncol(geb), ncol = length(mod_list))
rownames(coef_matrix_BIC) = colnames(geb)
colnames(coef_matrix_BIC) = mod_list



for (moi in mod_list) {
  
  print(moi)
  
  data = cbind(geb, mod[,moi])
  colnames(data)[ncol(data)] = "mark"
  data = as.data.frame(data)
  
  if (length(unique(data[,"mark"])) > 1) {
    
    # full model
    mod_tot = glm(mark ~ ., data = data, family = "binomial")
    # reduce full model using BIC
    mod_min = step(mod_tot, direction = "both", steps = 1e6, trace = F, k = log(nrow(data))) # k=log(num_obs) for BIC
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
    coef_matrix_BIC[model_min_coefs_sig_names,moi] = model_min_coefs[-1]
  }
  
}

# print coeficients in a heatmap
pdf("models_BIC.Me.pdf",height=20,width=20)
pheatmap(coef_matrix_BIC, color = cod.fun(31), breaks = seq(-10,10,length.out = 32), 
         cellwidth = 16, cellheight = 16, na_col = "grey",number_color = "darkgreen",
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = T,
         main="coefficients BIC regression, Me")
dev.off()


