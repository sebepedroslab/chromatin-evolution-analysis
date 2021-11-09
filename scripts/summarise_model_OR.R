glm_tables = function(model, null, ci=0.95, anova_test = "Chisq", model_name = "model", model_type = "glm") {
  
  # get model significance
  sig = anova(model, null, test = anova_test)
  
  # summarise model variables
  sum = summary(model)
  sut = as.data.frame(sum$coefficients)
  variables = rownames(sut)
  variables[1] = "Intercept"
  
  # get coefficients and odds ratio
  coefs = sut$Estimate
  
  # upper CI95 for ORs and coefficients
  coef_ciup = sut$Estimate + qnorm( 1 - (1-ci)/2 ) * sut$`Std. Error`
  coef_cido = sut$Estimate + qnorm(     (1-ci)/2 ) * sut$`Std. Error`
  
  # build output tables:
  # first, results from each variable 
  variables_table = data.frame(
    model = model_name,
    var = variables,
    coef = coefs,
    coef_CIdo = coef_cido,
    coef_CIup = coef_ciup,
    zscore = sut$`z value`,
    p = sut$`Pr(>|z|)`
  )
  
  # second, table summarising the entire model
  if (model_type == "glm") {
    num_obs = length(model$residuals)
    deviance = sum$deviance
    df_residual = sum$df.residual
    pval = sig[,"Pr(>Chi)"][2]
  } else if (model_type == "glmer") {
    num_obs = length(sum$residuals)
    deviance = sum$AICtab["deviance"]
    df_residual = sum$AICtab["df.resid"]
    pval = sig[,"Pr(>Chisq)"][2]
  } else {
    print("model_type must be either glm or glmer")
  }
  
  model_table = data.frame(
    model = model_name,
    num_obs = num_obs,
    anova_test = anova_test,
    p = pval,
    deviance_residual = deviance,
    df_residual = df_residual,
    CI = ci,
    variables = paste(variables[-1], collapse = ","),
    call = paste(capture.output(sum$call, file = NULL), collapse = "")
  )
  
  # output
  return(list("model_table" = model_table,
              "variable_table" = variables_table))
  
}

summarise_model_report_string = function(model) {
  # report model
  mod_sum = summary(model)
  mod_sum_out = as.data.frame(mod_sum$coefficients)
  # odds ratios:
  mod_sum_OR = 1/exp(mod_sum_out$Estimate)[-1] # odds ratio for each genotype (exclude intercept)
  mod_sum_ORcu = 1/exp(mod_sum$coefficients[-1,1] + qnorm(0.975) * mod_sum$coefficients[-1,2])
  mod_sum_ORcd = 1/exp(mod_sum$coefficients[-1,1] + qnorm(0.025) * mod_sum$coefficients[-1,2])
  mod_sum_var = rownames(mod_sum_out)[-1] # list of genotypes
  # create a string with odds ratios and CIs
  mod_sum_ORstring = c() ; van = 0
  for (vai in mod_sum_var) {
    van = van + 1
    mod_sum_ORstring[van] = sprintf("%s OR = %.3f ( %.3f - %.3f )", vai, signif(mod_sum_OR[van],5),  signif(mod_sum_ORcu[van],5), signif(mod_sum_ORcd[van],5) )
  }
  mod_sum_ORstring = paste(mod_sum_ORstring, collapse = "\n")
  return(mod_sum_ORstring)
}
