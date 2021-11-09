#### hygeofun ####

hygeofun = function(list_interest,annotation,gene_col,ano_col,outputname,
                    name_geneset,topnum,printfile=T,list_population=NA,
                    p_adj_method="BH",padjbool=T) {
  
  # Input 
  list_interest = unique(list_interest)
  annotation = annotation[,c(gene_col,ano_col)]
  colnames(annotation) = c("gene","annot")
  if (is.na(list_population)) {
    pol = as.vector(unlist(annotation$gene))
  } else {
    pol = unique(list_population)
  }
  
  # Set setting 
  ano_in_popl = subset(annotation, gene %in% pol)$annot
  ano_in_intl = subset(annotation, gene %in% list_interest & gene %in% pol)$annot
  
  if (length(ano_in_intl) > 0) {
    
    tab_in_popl = as.data.frame(table(ano_in_popl))
    tab_in_intl = as.data.frame(table(ano_in_intl))
    tab_in_all  = cbind(tab_in_popl,tab_in_intl[,2])
    tab_in_all  = tab_in_all[apply(tab_in_all!=0, 1, all),]
    colnames(tab_in_all) = c("annot","freq_in_pop","freq_in_int")
    tab_in_all$total_ano_in_pop = length(unlist(ano_in_popl))
    tab_in_all$total_ano_in_int = length(unlist(ano_in_intl))
    
    # Test 
    tab_in_all$pval = phyper(
      tab_in_all$freq_in_int-1,
      tab_in_all$freq_in_pop,
      tab_in_all$total_ano_in_pop-tab_in_all$freq_in_pop,
      tab_in_all$total_ano_in_int,
      lower.tail = F)
    
    if(padjbool) {
      tab_in_all$pval_adj = p.adjust(tab_in_all$pval,method = p_adj_method)
    } else {
      tab_in_all$pval_adj = tab_in_all$pval
    }
    tab_in_all          = tab_in_all[order(tab_in_all$pval_adj),]
    
    # Output 
    if(printfile){
      pdf(paste(outputname,".",name_geneset,".hygeo.",ano_col,".pdf",sep=""),height=4.5,width=4)
    }
    par(mar=c(5,12,5,2))
    ploti=barplot(height = rev(head(log(tab_in_all$pval_adj,10),topnum)),
                  names.arg = rev(head(tab_in_all$annot,topnum)),
                  xlim=c(0,-5),horiz=T,las=1,col="slategray3",border=NA,
                  cex.names=0.35,cex.axis=0.6,cex.lab=0.6,cex.sub=0.6,cex.main=0.6,
                  main=paste(name_geneset,ano_col),
                  sub =paste("n =",length(list_interest),
                             "genes | p adjust",p_adj_method,padjbool),
                  xlab="log(p)")
    abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
    abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
    text(x=0,ploti,labels = paste("p =",signif(rev(head(tab_in_all$pval_adj,topnum)),3),"n =",rev(head(tab_in_all$freq_in_int,topnum))),
         col="red",pos=4,cex=0.35)
    
    if(printfile){
      dev.off()
      write.table(tab_in_all,
                  file=paste(outputname,".",name_geneset,".hygeo.",ano_col,".txt",sep=""),
                  row.names = F,sep="\t",quote = F)
    }
  } else {
    print("skip, no annotations in interest list!")
  }
  
}




#### TOPGO ####

topgofun  = function(list_interest,gomap,outputname,name_geneset,ontologyset,tg_test,tg_algorithm,topnum=20,nodesize=10,printfile=T,p_adj_method="BH") {
  
  library(topGO)
  
  # Input 
  list_interest = unique(list_interest)
  genom = names(gomap)
  gesel = factor(as.integer(genom %in% list_interest))
  names(gesel) = genom
  
  # shortened go mappings without empty transcripts
  gomap_nonempty = gomap[lapply(gomap,length)>0]
  
  if(printfile){
    pdf(file=paste(outputname,".",name_geneset,".topgo",".",tg_test,tg_algorithm,".pdf",sep=""),height=4.5,width=4)
  }
  par(mar=c(5,12,5,2))
  
  topgo_tau_tot = data.frame()
  
  if (length(list_interest[list_interest %in% names(gomap_nonempty)])>1) {
    
    for (ontologyseti in ontologyset) {
      # topGO setup 
      
      GOdata = new("topGOdata", ontology=ontologyseti, allGenes=gesel,
                   annot=annFUN.gene2GO, gene2GO=gomap)
      
      num_interest_feasible = sum(GOdata@feasible & genom %in% list_interest)
      
      # topGO anàlisi 
      topgo_res = runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test)
      topgo_tau = GenTable(GOdata, pval_test = topgo_res,
                           orderBy = "pval_test", 
                           topNodes = length(usedGO(object = GOdata)))
      topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
      topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
      topgo_tau$ontology = ontologyseti
      topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
      
      # Output 
      ploti=barplot(height = rev(head(log(topgo_tau$pval_test,10),topnum)),
                    names.arg = rev(head(paste(topgo_tau$Term,topgo_tau$GO.ID),topnum)),
                    xlim=c(0,-5),horiz=T,las=1,col="slategray3",border=NA,
                    cex.names=0.35,cex.axis=0.6,cex.lab=0.6,cex.sub=0.6,cex.main=0.6,
                    main=paste(name_geneset,"top GO:",ontologyseti,tg_test,tg_algorithm),
                    sub =paste("n=",num_interest_feasible,"/",length(list_interest), sep=""),
                    xlab="log(p)")
      abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
      abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
      text(x=0,ploti,labels = paste("p =",signif(rev(head(topgo_tau$pval_test,topnum)),3)),
           col="red",pos=4,cex=0.35)
    }
    
  }else {
    print("skip, no annotations in interest list!")
  }
  
  if(printfile){
    write.table(
      topgo_tau_tot,
      file=paste(outputname,".",name_geneset,".topgo",".",tg_test,tg_algorithm,".txt",sep=""),
      sep="\t", quote=F, col.names=T, row.names=F, append = F)
    dev.off()
  }
  
}



#### VENN DIAGRAMS ####

venn.two = function(list1,list2,catname1,catname2,main="Venn",
                    col1="green3",col2="magenta3",
                    eulerbool=T,print.mode=c("raw","percent")) {
  
  library(VennDiagram)
  library(grid)
  
  # compute set overlaps, intersections, etc
  list_uniq      = base::unique(c(list1, list2))
  list_intersect = base::intersect(list1, list2)
  list_diff_i    = base::setdiff(list1, list2)
  list_diff_j    = base::setdiff(list2, list1)
  
  # draw venn
  venn = draw.pairwise.venn(
    area1 = length(list1),
    area2 = length(list2),
    cross.area = length(list_intersect),
    category=c(catname1,catname2),
    fontfamily="Helvetica",
    cat.fontfamily = "Helvetica",
    col=c(col1,col2),ext.text = F,
    cat.col=c(col1,col2),
    ind = F, print.mode = print.mode,
    euler.d = eulerbool)
  grid::grid.newpage()
  grid::grid.draw(venn)
  grid::grid.text(main,x=0.5,y=0.95,
                  gp = gpar(fontsize = 8, fontface = "bold"))
  
  # output
  output = list(
    "catname_i"      = catname1,
    "catname_j"      = catname2,
    "title"          = main,
    "list_uniq"      = list_uniq,
    "list_intersect" = list_intersect,
    "list_diff_i"    = list_diff_i,
    "list_diff_j"    = list_diff_j,
    "venn" = venn
    )
  
  return(output)
  
}


venn.three = function(list1,list2,list3,
                      catname1,catname2,catname3,
                      main="Venn",
                      col1="green3",col2="magenta3",col3="orange",
                      eulerbool=T,print.mode=c("raw","percent")) {
  
  library(VennDiagram)
  library(grid)
  
  # compute set overlaps, intersections, etc
  intersect_n12  = base::intersect(list1, list2)
  intersect_n13  = base::intersect(list1, list3)
  intersect_n23  = base::intersect(list2, list3)
  intersect_n123 = base::intersect(intersect_n12, intersect_n13)
  
  # draw venn
  venn = draw.triple.venn(
    area1 = length(list1),
    area2 = length(list2),
    area3 = length(list3),
    n12 = length(intersect_n12),
    n13 = length(intersect_n13),
    n23 = length(intersect_n23),
    n123 = length(intersect_n123),
    category=c(catname1,catname2,catname3),
    fontfamily="Helvetica",
    cat.fontfamily = "Helvetica",
    col=c(col1,col2,col3),ext.text = F,
    cat.col=c(col1,col2,col3),
    ind = F, print.mode = print.mode,
    euler.d = eulerbool, scaled = T)
  grid::grid.newpage()
  grid::grid.draw(venn)
  grid::grid.text(main,x=0.5,y=0.95,
                  gp = gpar(fontsize = 8, fontface = "bold"))
  
  # output
  output = list(
    "venn" = venn,
    "catname_1"      = catname1,
    "catname_2"      = catname2,
    "catname_3"      = catname3,
    "intersect_12" = intersect_n12,
    "intersect_13" = intersect_n13,
    "intersect_23" = intersect_n23,
    "intersect_123" = intersect_n123
  )
  
  return(output)
  
}


