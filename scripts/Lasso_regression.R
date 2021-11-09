lasso_regression_marks=function(variable_matrix,predictor_matrix){  #reduce predictor matrix if necessary to shared species!
	
  library(glmnet)
  
	variable_matrix[is.na(variable_matrix)]=0  #remove NAs
	
	a=names(which(colSums(variable_matrix)<(nrow(variable_matrix)-1))) ##we cannot regress  if a mark is present everywhere or everywhere minus one.
	b=names(which(colSums(variable_matrix)>1)) #we cannot regress if a mark present only in 1 sps.
	to_keep=intersect(a,b)
	
	variable_matrix=variable_matrix[,to_keep]
	
	
	vars=colnames(variable_matrix)
	out.m=matrix(ncol=ncol(variable_matrix),nrow=ncol(predictor_matrix))
	colnames(out.m)=colnames(variable_matrix)
	rownames(out.m)=colnames(predictor_matrix)
	
	for(v in vars){
		#browser()
		m=cbind.data.frame(variable_matrix[,v],predictor_matrix)
		colnames(m)=c(v,colnames(predictor_matrix))
		x=as.matrix(m[,-1])
		y=as.vector(m[,v])
		cv.lasso=cv.glmnet(x,y,alpha=1,family="binomial")
		model=glmnet(x,y,alpha=1,family="binomial",lambda=cv.lasso$lambda.min)
		
		out.m[,v]=coef(model)[rownames(out.m),]
	}
	
	return(out.m)
}



lasso_regression_single_mark=function(mark,variable_matrix,predictor_matrix){
	variable_matrix[is.na(variable_matrix)]=0 
	
	m=cbind.data.frame(variable_matrix[,mark],predictor_matrix)
	colnames(m)=c(mark,colnames(predictor_matrix))
	
	x=as.matrix(m[,-1])
	y=as.vector(m[,mark])
	cv.lasso=cv.glmnet(x,y,alpha=1,family="binomial")
	model=glmnet(x,y,alpha=1,family="binomial",lambda=cv.lasso$lambda.min)
	
	return(model)
}