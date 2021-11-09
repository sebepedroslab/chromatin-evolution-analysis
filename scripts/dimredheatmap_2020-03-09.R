dimredfun = function(
  matriu, # input data (can be gene expression, presence/absence, etc); must be numeric matrix with rownames and colnames
  outputname = "output", # prefix for the output file
  varname = "title", # plot title
  cols_are = "samples", # title for cols 
  rows_are = "genes", # title for rows
  isbidi = F, # use a bidirectional color map (red-to-white-to-bue) or a continuous color map (white-to-blue)?
  cols_dist_method = "spearman", # dist method for column clustering ; works with distance metrics ("euclidean" in dist() function) or correlation metrics ("spearman", "kendall", "pearson" in cor() function).
  rows_dist_method = "spearman", # id for rows
  clus_method = "ward.D2",  # clustering method for hclust; e.g. "ward.D2"
  printpdfpca = T, # print PCAs?
  printpdfheatmaps = T,  # print heatmaps?
  hm_height = 8,  # heatmap height (pdf)
  hm_width = 8, # heatmap width (pdf)
  omit_clustering_if_n_observations = 1000 # omit clustering in rows and cols if there are more than X observations
  ) {
  
  
  mi = matriu
  
  library(gplots)
  library(ape)
  # library(factoextra)
  
  graphics.off()
  
  # Colors
  if (isbidi == T) {
    hop = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
    div = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
    
  } else {
    hop = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
    div = colorRampPalette(c("firebrick4","orangered",
                             "floralwhite",
                             "deepskyblue","dodgerblue4"))
  }
  
  # PCA
  message(paste("# PCA"))
  pic                  = prcomp(t(mi))
  pic$variancefraction = pic$sdev^2/sum(pic$sdev^2)
  
  # Correlation, distance & clustering
  message(paste("# Cols",cols_dist_method,clus_method))
  if (cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") {
    mic.cor  = cor(mi,method = cols_dist_method)
    mic.dist = as.dist(1-mic.cor)
    mic.clus = hclust(mic.dist,method=clus_method)
  } else {
    mic.cor  = data.frame()
    mic.dist = dist(t(mi),method=cols_dist_method)
    mic.clus = hclust(mic.dist,method=clus_method)
  }
  
  message(paste("# Rows",rows_dist_method,clus_method))
  if (rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") {
    mir.cor  = cor(t(mi),method=rows_dist_method)
    mir.dist = as.dist(1-mir.cor)
    mir.clus = hclust(mir.dist,method=clus_method)
  } else {
    mir.cor  = data.frame()
    mir.dist = dist(mi,method=rows_dist_method)
    mir.clus = hclust(mir.dist,method=clus_method)
  }
  
  # PCoA
  message(paste("# PCoA"))
  pio      = pcoa(mic.dist)
  pio12    = as.data.frame(pio$vectors[,c(1,2,3)])
  pio12$sp = rownames(pio12)
  
  # Plots PCA/PCOA
  if (printpdfpca) {
    
    message(paste("# Plot PCA & PCoA"))
    pdf(paste(outputname,".reddim.pdf",sep=""),height=9,width=8)
    par(mfrow=c(2,2))
    
    # Plot PCA: 12
    plot(pic$x[,c(1,2)],col="red",main="PCA 1&2",
         xlab=paste("PCA 1",signif(pic$variancefraction[1]*100,digits=3),"% variance"),
         ylab=paste("PCA 2",signif(pic$variancefraction[2]*100,digits=3),"% variance"))
    text(pic$x[,c(1,2)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: 13
    plot(pic$x[,c(1,3)],col="red",main="PCA 1&3",
         xlab=paste("PCA 1",signif(pic$variancefraction[1]*100,digits=3),"% variance"),
         ylab=paste("PCA 3",signif(pic$variancefraction[3]*100,digits=3),"% variance"))
    text(pic$x[,c(1,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: 23
    plot(pic$x[,c(2,3)],col="red",main="PCA 2&3",
         xlab=paste("PCA 2",signif(pic$variancefraction[2]*100,digits=3),"% variance"),
         ylab=paste("PCA 3",signif(pic$variancefraction[3]*100,digits=3),"% variance"))
    text(pic$x[,c(2,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCA: eigenvalues
    plot(pic$variancefraction,type="b",col="red",main="PCA variance")
    text(pic$variancefraction,labels=signif(pic$variancefrac,digits=3),pos=3,col="black")
    
    # Plot PCoA: 12
    plot(x=pio12$Axis.1,y=pio12$Axis.2,col="blue",main="PCoA 1&2",
         xlab=paste("Axis 1",signif(pio$values$Relative_eig[1]*100,digits=3),"% variance"),
         ylab=paste("Axis 2",signif(pio$values$Relative_eig[2]*100,digits=3),"% variance"))
    text(pio$vectors[,c(1,2)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: 13
    plot(x=pio12$Axis.1,y=pio12$Axis.3,col="blue",main="PCoA 1&3",
         xlab=paste("Axis 1",signif(pio$values$Relative_eig[1]*100,digits=3),"% variance"),
         ylab=paste("Axis 3",signif(pio$values$Relative_eig[3]*100,digits=3),"% variance"))
    text(pio$vectors[,c(1,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: 23
    plot(x=pio12$Axis.2,y=pio12$Axis.3,col="blue",main="PCoA 2&3",
         xlab=paste("Axis 2",signif(pio$values$Relative_eig[2]*100,digits=3),"% variance"),
         ylab=paste("Axis 3",signif(pio$values$Relative_eig[3]*100,digits=3),"% variance"))
    text(pio$vectors[,c(2,3)],labels=factor(colnames(mi)),pos=3,col="black")
    # Plot PCoA: eigenvalues
    plot(pio$values$Relative_eig,col="blue",main="PCoA relative eigenvalues",type="b",xlab="Axis rank",ylab="Fraction variance explained")
    text(pio$values$Relative_eig,labels=signif(pio$values$Relative_eig,digits=3),pos=3,col="black")
    lines(pio$values$Broken_stick[pio$values$Broken_stick>.01],
          col="gray",main="PCoA variance broken stick",type="b")
    
    # Plot clustering of columns
    plot(mic.clus, main = paste(cols_are," (",ncol(mi),")",sep=""), sub = paste("clustering:",clus_method,"distances:",cols_dist_method))
    
    dev.off()
  }
  
  
  
  # plot heatmaps
  message(paste("# Plot heatmaps"))
  if (printpdfheatmaps){
    
    pdf(paste(outputname,".heatmaps.pdf",sep=""),height=hm_height,width=hm_width)
    par(mfrow=c(1,1))
    
    labrowbool = if(nrow(mi) > 500) { F } else { rownames(mi) }
    labcolbool = if(ncol(mi) > 500) { F } else { colnames(mi) }
    
    # heatmap samples vs variables
    heatmap.2(mi,
              col =hop(31),dendrogram = "both",symkey=F,trace="none",scale="none",
              useRaster=F,keysize = 1,main = paste(cols_are," (",ncol(mi),")"," ~ ",
                                                   rows_are," (",nrow(mi),")",sep=""),symm=F,
              labRow = labrowbool,
              labCol = labcolbool,
              Colv=as.dendrogram(mic.clus),
              Rowv=as.dendrogram(mir.clus),
              margins=c(10,10),na.color = "grey90",key.title = varname)
    
    # heatmap variables correlation
    if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < omit_clustering_if_n_observations ) {
      
      labrowbool = if(nrow(mir.cor) > 500) { F } else { rownames(mir.cor) }
      labcolbool = if(ncol(mir.cor) > 500) { F } else { colnames(mir.cor) }
      
      heatmap.2(mir.cor,
                col =div(31),dendrogram = "both",symbreaks=T,trace="none",scale="none",
                useRaster=F,keysize = 1,main = paste("Pairwise",rows_are),symm=F,
                labRow = labrowbool,
                labCol = labcolbool,
                Colv=as.dendrogram(mir.clus),
                Rowv=as.dendrogram(mir.clus),
                margins=c(10,10),na.color = "grey90",key.title = "Corr. coef.")
    }
    
    # heatmap samples correlation
    if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < omit_clustering_if_n_observations ) {
      
      labrowbool = if(nrow(mic.cor) > 500) { F } else { rownames(mic.cor) }
      labcolbool = if(ncol(mic.cor) > 500) { F } else { colnames(mic.cor) }
      
      heatmap.2(mic.cor,
                col =div(31),
                dendrogram = "both",symbreaks=T,trace="none",scale="none",
                useRaster=F,keysize = 1,main = paste("Pairwise",cols_are),symm=F,
                Colv=as.dendrogram(mic.clus),
                labRow = labrowbool,
                labCol = labcolbool,
                Rowv=as.dendrogram(mic.clus),
                margins=c(10,10),na.color = "grey90",key.title = "Corr. coef.")
    }
    
    dev.off()
  }
  
  # output: clusterings etc.
  return.list = list(
    "pio"      = pio,
    "pic"      = pic,
    "mic.clus" = if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < omit_clustering_if_n_observations ) { mic.clus } else { NA },
    "mic.dist" = if ((cols_dist_method == "spearman" | cols_dist_method == "pearson" | cols_dist_method == "kendall") & nrow(mic.cor) < omit_clustering_if_n_observations ) { mic.dist } else { NA },
    "mir.clus" = if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < omit_clustering_if_n_observations ) { mir.clus } else { NA },
    "mir.dist" = if ((rows_dist_method == "spearman" | rows_dist_method == "pearson" | rows_dist_method == "kendall") & nrow(mir.cor) < omit_clustering_if_n_observations ) { mir.dist } else { NA }
  )
  
  return(return.list)
  
}
