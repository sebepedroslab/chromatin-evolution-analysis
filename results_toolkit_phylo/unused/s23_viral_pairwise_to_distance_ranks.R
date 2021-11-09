library(ape)
library(phytools)
library(scales)

domlist=c("SET","BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","HIF-1","Histone",
          "LinkerHistone","PHD","PTIP","SIR2","SNF2_N","TUDOR",
          "zf-CCHH","zf-CXXC","Hist_deacetyl","Acetyltransf_1")

# domlist=c("DOT1")


for (dom in domlist) {
  
  print(dom)
  
  # define input
  ali_fn = sprintf("results_viruses/pairwise_alignments/vir.%s.diamond.to_cellular.diamond.csv.gz", dom)
  tre_fn = sprintf("results_viruses/pairwise_alignments/vir.%s.diamond.newick", dom)
  
  # load
  ali = read.table(
    gzfile(ali_fn), header = F, 
    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen",
                  "qstart","qend","sstart","send","evalue","bitscore"))
  
  # filter out secondary alignments of each hit
  ali = ali [!duplicated(ali[c("qseqid","sseqid")]),]
  
  # pairwise bitscore matrices
  alm_score = reshape(ali[,c("qseqid","sseqid","bitscore")], idvar = "qseqid", timevar = "sseqid", direction = "wide")
  ali_sequences = alm_score$qseqid
  alm_score = as.matrix(alm_score[,-1])
  rownames(alm_score) = ali_sequences
  colnames(alm_score) = ali_sequences
  
  # if bitscore is NA (no alignment for that pair), change to minimum observed value
  alm_score_min_val = apply(alm_score, 1, function(x) min(x, na.rm = T))
  alm_score_ixs_noa = apply(alm_score, 1, function(x) which(is.na(x)))
  alm_score_ixs_ali = apply(alm_score, 1, function(x) which(!is.na(x)))
  
  for (i in seq(nrow(alm_score))) {
    alm_score[i,alm_score_ixs_noa[[i]]] = alm_score_min_val[i]
  }
  
  # normalised bitscore: value divided by max value in the dataset (including self?)
  alm_score = alm_score / apply(alm_score, 1, max)
  
  # hit-weighted bitscore: account for the number of alignments for the query and query+hit
  # hits to query
  alm_hits_qu = lapply(alm_score_ixs_ali, function(x) length(names(x)))
  # hits to both subject and query
  alm_hits_qh_intersection = lapply(alm_score_ixs_ali, function(x) 
    lapply(x, function(z) length(intersect(names(x), names(alm_score_ixs_ali[[z]]))) )
  )
  # pairwise matrix
  alm_qhint = matrix(0, nrow = nrow(alm_score), ncol = ncol(alm_score))
  rownames(alm_qhint) = rownames(alm_score)
  colnames(alm_qhint) = colnames(alm_score)
  for (i in seq(nrow(alm_qhint))) {
    alm_qhint[i,names(alm_hits_qh_intersection[[i]])] = as.numeric(alm_hits_qh_intersection[[i]]) / alm_hits_qu[[i]]
  }
  
  # weight scores
  alm_score = alm_score * alm_qhint
  
  
  # find indexes of viruses
  ix_vir = grep(rownames(alm_score), pattern = "^vir_")
  ix_cel = grep(rownames(alm_score), pattern = "^vir_", invert = T)
  
  if (length(ix_vir)>1){
    pdf(sprintf("results_viruses/pairwise_distributions-%s.pdf", dom), width = 6, height = 4)
    # scores from viruses to cells or viruses
    alm_score_vir2cel = alm_score[ix_vir,ix_cel]
    alm_score_vir2vir = alm_score[ix_vir,ix_vir]
    
    # closest cell to each vir
    hit_score_vir2cel_ixs = apply(alm_score_vir2cel, 1, which.max)
    hit_score_vir2cel_val = apply(alm_score_vir2cel, 1, max)
    # closest vir to each vir
    diag(alm_score_vir2vir) = -Inf
    hit_score_vir2vir_ixs = apply(alm_score_vir2vir, 1, which.max)
    hit_score_vir2vir_val = apply(alm_score_vir2vir, 1, max)
    
    hit_vir = data.frame(
      gene = rownames(alm_score)[ix_vir],
      hit_vir_score = hit_score_vir2vir_val,
      hit_vir_gene  = colnames(alm_score_vir2vir)[hit_score_vir2vir_ixs],
      hit_cel_score = hit_score_vir2cel_val,
      hit_cel_gene  = colnames(alm_score_vir2cel)[hit_score_vir2cel_ixs]
    )
    hit_vir$l2r_scores = log2(hit_vir$hit_vir_score / hit_vir$hit_cel_score)
    # plot(sort(hit_vir$l2r_scores), col="blue")
    # abline(h=0, lty=2)
    
    
    
    # plot distribution of hits for each non-euk best hit
    hit_vir_cels = hit_vir$hit_cel_gene
    cel_hit_done = c()
    # cel_hit_list = list(euk=c(),vir=c(),arc=c(),bac=c())
    for (i in seq(length(hit_vir_cels))) {
      cel_hit = hit_vir[i,"hit_cel_gene"]
      cel_hit_tax = strsplit(cel_hit, split = "_")[[1]][1]
      if ( !cel_hit_tax %in% c("arc","bac","vir")) {  cel_hit_tax = "euk" }
      if (!cel_hit %in% cel_hit_done) {
        vir_hit = hit_vir[i,"hit_vir_gene"]
        vir_hit_score  = hit_vir[i,"hit_vir_score"]
        cel_hit_scores = sort(alm_score[,cel_hit], decreasing = T)
        cel_hit_scores_col = rep("slategray", length(cel_hit_scores))
        cel_hit_scores_col[grepl(names(cel_hit_scores), pattern = "^vir_")] = "deeppink"
        cel_hit_scores_col[grepl(names(cel_hit_scores), pattern = "^bac_")] = "blue"
        cel_hit_scores_col[grepl(names(cel_hit_scores), pattern = "^arc_")] = "darkturquoise"
        cel_hit_scores_tax = rep("euk", length(cel_hit_scores))
        cel_hit_scores_tax[grepl(names(cel_hit_scores), pattern = "^vir_")] = "vir"
        cel_hit_scores_tax[grepl(names(cel_hit_scores), pattern = "^bac_")] = "bac"
        cel_hit_scores_tax[grepl(names(cel_hit_scores), pattern = "^arc_")] = "arc"
        cel_hit_scores_tax_order = cel_hit_scores_tax[!duplicated(cel_hit_scores_tax)]
        cel_hit_scores_tax_vrank = which(cel_hit_scores_tax_order == "vir")
        plot(cel_hit_scores, col=cel_hit_scores_col, 
             main=sprintf("viral rank = %i | %s",cel_hit_scores_tax_vrank,dom), 
             sub = sprintf("%s (hit from %s)", cel_hit, vir_hit),
             ylab="alignment score", cex.sub=0.8)
        abline(h=vir_hit_score, lty=2, col="orange")
        legend("topright",legend=c("euk","vir","arc","bac"), col = c("slategray","deeppink","darkturquoise","blue"), pch=1)
        cel_hit_done = c(cel_hit_done,cel_hit)
        
        # # store values
        # cel_hit_list[["euk"]] = c(cel_hit_list[["euk"]], cel_hit_scores[cel_hit_scores_tax == "slategray" & cel_hit_scores > 0])
        # cel_hit_list[["vir"]] = c(cel_hit_list[["vir"]], cel_hit_scores[grepl(names(cel_hit_scores), pattern = "^vir_") & cel_hit_scores > 0])
        # cel_hit_list[["arc"]] = c(cel_hit_list[["arc"]], cel_hit_scores[grepl(names(cel_hit_scores), pattern = "^arc_") & cel_hit_scores > 0])
        # cel_hit_list[["bac"]] = c(cel_hit_list[["bac"]], cel_hit_scores[grepl(names(cel_hit_scores), pattern = "^bac_") & cel_hit_scores > 0])
        
      }
    }
    
    # boxplot(cel_hit_list, col=c("slategray","deeppink","darkturquoise","blue"),
    #         main="global",outline = F, ylab="Distance")
    # plot(NA,NA,xlim = c(0,1), ylim=c(0,10), col="slategray", xlab="Distance", ylab="Density")
    # if (length(cel_hit_list[["euk"]]) > 1) {lines(density(cel_hit_list[["euk"]]), col="slategray")}
    # if (length(cel_hit_list[["vir"]]) > 1) {lines(density(cel_hit_list[["vir"]]), col="deeppink")}
    # if (length(cel_hit_list[["arc"]]) > 1) {lines(density(cel_hit_list[["arc"]]), col="darkturquoise")}
    # if (length(cel_hit_list[["bac"]]) > 1) {lines(density(cel_hit_list[["bac"]]), col="blue")}
    
    dev.off()
  }
  
  
  #### PLOT TREE ####
  
  # many cells will become zero because they have zero overlap
  # if this happens, add a min value derived from the alignment length
  # alm_score_weighted[alm_score_weighted==0] = mean(alm_score_weighted) / max(ali$length)
  # score to distance
  dis = -log(alm_score)
  dis[is.infinite(dis)] = max(dis[!is.infinite(dis)])+1
  dis = as.dist(dis)
  
  
  # # plot tree
  # pdf(sprintf("%s.pdf", tre_fn), width = 12, height = 12)
  # 
  # # save as tree  
  # phyl = ape::as.phylo(hclust(dis, method = "ward.D2"))
  # ape::write.tree(phyl, file=tre_fn)
  # 
  # phyl = midpoint.root(phyl)
  # phyl = ladderize(phyl)
  # 
  # # dataframe of edges
  # phyl_edge           = as.data.frame(phyl$edge)
  # colnames(phyl_edge) = c("edge_start","edge_end")
  # phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
  # phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
  # # dataframe of nodes
  # phyl_nods = data.frame(node = c(phyl$tip.label, phyl$node.label))
  # phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
  # phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)
  # 
  # # taxa in each node
  # phyl_nods$taxa = ""
  # phyl_nods[phyl_nods$is_tip,"taxa"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "_", simplify = T)[,1]
  # # phyl_nods[!phyl_nods$taxa %in% c("arc","bac","vir"),"taxa"] = "euk"
  # # gene name
  # phyl_nods$gene = ""
  # phyl_nods[phyl_nods$is_tip,"gene"] = stringr::str_split(phyl_nods[phyl_nods$is_tip,"node"],pattern = "\\|", simplify = T)[,1]
  # 
  # # merge them
  # phyl_edge = merge(phyl_edge, phyl_nods, all.x = T, all.y = F, by.x = "edge_end", by.y = "edge_end",sort = F)
  # 
  # # taxa in each edge
  # phyl_edge$taxa = stringr::str_split(phyl_edge$node,pattern = "_", simplify = T)[,1]
  # # phyl_edge[!phyl_edge$taxa %in% c("arc","bac","vir"),"taxa"] = "euk"
  # 
  # # add colors
  # phyl_edge$color = "slategray3"
  # phyl_edge[phyl_edge$taxa %in% c("bac"), "color"] = "blue"
  # phyl_edge[phyl_edge$taxa %in% c("arc"), "color"] = "darkturquoise"
  # phyl_edge[phyl_edge$taxa %in% c("vir"), "color"] = "deeppink"
  # 
  # # add colors
  # phyl_nods$color = "slategray4"
  # phyl_nods[phyl_nods$taxa %in% c("bac"), "color"] = "blue"
  # phyl_nods[phyl_nods$taxa %in% c("arc"), "color"] = "darkturquoise"
  # phyl_nods[phyl_nods$taxa %in% c("vir"), "color"] = "deeppink"
  # 
  # # plot unrooted
  # plot.phylo(phyl, font=1, type="u", edge.color = "slategray3", root.edge = T, show.tip.label = F, underscore = T)
  # tiplabels(col = phyl_nods$color,frame = "none", pch = 19, cex=0.8)
  # tiplabels(text = phyl_nods[phyl_nods$taxa=="vir","gene"],tip = which(phyl_nods$taxa == "vir"), frame="none", col=alpha("deeppink2",0.6), cex=0.4)
  # add.scale.bar()
  # title(sprintf("%s\nn=%i sequences",dom,length(phyl$tip.label)))
  # 
  # legend("topright", col=c("slategray4","blue","darkturquoise","deeppink"), pch=19, legend = c("eukaryotes","bacteria","archaea","viruses"))
  # 
  # dev.off()  
}

