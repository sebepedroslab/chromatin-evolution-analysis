# libraries
library("stringr")
library("ape")
library("readxl")
library("vioplot")
library("scales")

# input
ref_fn = "human_histones.dict_classes.csv"
pro_fn = "all.Histone.to_proteomics.csv.gz" # histones that appear in the proteomics datasets
mod_fn = "../../results_PTMs/consensus_modifications_perseq.xlsx" # list of canonical histones used in proteomics analyse
ali_fn = c("euk.Histone.domains.diamond_to_all.csv.gz","euk.Histone.domains.diamond_to_human.csv.gz","euk.Histone.domains.diamond_to_histdb.csv.gz","vir.Histone.domains.diamond_to_histdb.csv.gz","arc.Histone.domains.diamond_to_histdb.csv.gz","all.Histone.domains.diamond_to_histdb.csv.gz")
ali_fn = c("all.Histone.to_histdb.csv.gz","euk.Histone.to_histdb.csv.gz","euk.Histone.to_human.csv.gz")
his_list = c("H2A","H2B","H3","H4","macroH2A","H2AZ","other_CENP")

# ali_fn = c("euk.Histone.domains.diamond_to_histdb.csv")
# his_list = c("H2AZ")

ref = read.table(ref_fn, sep = "\t", col.names = c("gene","class"))
ref$class = factor(ref$class)
class_col = data.frame(histone = levels(ref$class))
class_col$color = rainbow(n=nrow(class_col), v=0.9, end = 0.8)

# load list of histones that appear in the proteomics dataset
pro = read.table(pro_fn, col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
# keep only identical hits
pro = pro [ pro$pident == 100, ]

# list of canonical histones identified in the proteomics dataset
# keep only reference sequences that appear in the proteomics dataset
his_in_pro = list()
hip_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")
hip_cols = c("blue4","darkorange3","darkgreen","darkmagenta","darkolivegreen3","darkolivegreen4")
for (hip in hip_list) {
     
     mod = read_excel(mod_fn, sheet=hip)
     mod = mod [ !is.na(mod$Positions_in_Reference), ]
     his_in_mod = gsub("\\|.*","",mod$Protein_ID)
     his_in_pro[[hip]] = pro [ pro$sseqid %in% his_in_mod , "qseqid" ]
     
}

for (ali_fi in ali_fn) {
    
    print(ali_fi)
    
    ali_fi = stringr::str_split(ali_fi, pattern = ".csv.gz")[[1]][1]
    
    # load
    ali = read.table(
        gzfile(sprintf("%s.csv.gz",ali_fi)), sep = "\t", 
        col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
        stringsAsFactors = FALSE)
    
    # retrieve matched histone
    ali$histone = as.factor(stringr::str_split(ali$sseqid, pattern = "\\|",simplify = TRUE)[,2])
    ali = merge(ali, class_col, all.x=TRUE, all.y = FALSE, by.x="histone", by.y = "histone")
    ali = ali[order(ali$qseqid,ali$histone,ali$bitscore,decreasing=TRUE),]
    
    # pdf reports
    pdf(sprintf("%s.diagnostic.pdf",ali_fi), width = 20, height = 8)
    layout(mat = matrix(c(1:10),nrow = 2))
    for (his in his_list) {
        
        # retrieve hits to this histone
        aii = ali[ali$histone == his,]
        aii_main_hits = aii[,c("qseqid","histone","pident","evalue","bitscore","color")]
        aii_main_hits = aii_main_hits[order(aii_main_hits$bitscore,decreasing = TRUE),]
        aii_main_hits = aii_main_hits[!duplicated(aii_main_hits[,c("qseqid","histone")]),]
        
        # plot main hits pident dist
        hist(aii$pident, breaks = 40 ,col = "gray85", border = "white", main = his, xlim = c(0,100), xlab = "pident", ylim = c(0,3600))
        # add dots with values of hits to proteomics histones
        m=0
        for (n in 1:length(hip_list)) {
             hip = hip_list[n]
             hic = hip_cols[n]
             aii_hip = aii_main_hits [ aii_main_hits$qseqid %in% his_in_pro[[hip]] , "pident" ]
             if (length(aii_hip) > 0 ) {
               # points(x=aii_hip, y = rnorm(length(aii_hip), mean = n * 100, sd = 10), col = alpha(hic,0.6))
               vioplot(aii_hip, at = 3000 + m * 250, horizontal = TRUE, col = hic, ylim = c(0,100), add = TRUE, wex = 200)
               m = m + 1
             }
        }
        legend("topleft",hip_list, col = hip_cols, pch = 1, bty = "n")
        
     #    # boxplots instead of points (wider)
     #    plot(NA,NA, xlim = c(0,100), ylim = c(0,10), xlab = "pident")
     #    for (n in 1:length(hip_list)) {
     #         hip = hip_list[n]
     #         hic = hip_cols[n]
     #         aii_hip = aii_main_hits [ aii_main_hits$qseqid %in% his_in_pro[[hip]] , "pident" ]
     #         if (length(aii_hip) > 0 ) {
     #           vioplot(aii_hip, at = n, horizontal = TRUE, col = hic, ylim = c(0,100), add = TRUE)
     #         }
     #    }
     #    legend("topleft",hip_list, col = hip_cols, pch = 1, bty = "n")
        
        # among the hits to this histone, find which other histones have the hit
        aii_queries = unique(aii$qseqid)
        aii_other_hits = ali[ali$qseqid %in% aii_queries & ali$histone != his,c("qseqid","histone","pident","evalue","bitscore","color")]
        aii_other_hits = aii_other_hits[order(aii_other_hits$bitscore,decreasing = TRUE),]
        aii_other_hits = aii_other_hits[!duplicated(aii_other_hits[,c("qseqid","histone")]),]
        
        aii_main_table = table(aii_main_hits$histone) + table(aii_other_hits$histone)
        pie(aii_main_table, labels = paste(names(aii_main_table),"\nn =",aii_main_table), 
            main=sprintf("Num alignments = %i", sum(aii_main_table)),
            cex=0.8)
        
        aii_other_table = table(aii_other_hits$histone)
        if (sum(aii_other_table) > 0) {
            
            # merge
            aii_merged = merge(aii_main_hits, aii_other_hits, by="qseqid", all.y = TRUE, all.x = TRUE)
            aii_merged$histone.y = droplevels(aii_merged$histone.y)
            aii_merged$pident.xy = aii_merged$pident.x / aii_merged$pident.y
            aii_merged$bitscore.xy = aii_merged$bitscore.x / aii_merged$bitscore.y
            aii_merged$evalue.xy = aii_merged$evalue.x / aii_merged$evalue.y
            
            # plot pident    
            plot(aii_merged$pident.x, aii_merged$pident.y, 
                 col=aii_merged$color.y,
                 xlim=c(0,100), ylim=c(0,100),
                 xlab = sprintf("pident to %s", his),
                 ylab = "pident to another histone alignment",
                 main = sprintf("pident %s",his))
            abline(a=0, b=1, lty=2, col="gray")
            legend("bottomleft", 
                   legend = levels(aii_merged$histone.y), 
                   col = class_col[as.character(class_col$histone) %in% levels(aii_merged$histone.y),"color"], 
                   bty="n", cex=0.8, pch=1, )
            # ratio dist
            hist(aii_merged$bitscore.xy, breaks = 40, xlim = c(0,5),
                 col = "gray85", border = "white", main = his,xlab = sprintf("pident %s / pident secondary", his)) 
            abline(v=1, col="gray", lty=2)
            
            # plot bitscore
            plot(aii_merged$bitscore.x, aii_merged$bitscore.y,
                 col=aii_merged$color.y,
                 xlim=c(0,300), ylim=c(0,300),
                 xlab = sprintf("bitscore to %s", his),
                 ylab = "bitscore to another histone alignment",
                 main = sprintf("bitscore %s",his))
            abline(a=0, b=1, lty=2, col="gray")
            legend("bottomright", 
                   legend = levels(aii_merged$histone.y), 
                   col = class_col[as.character(class_col$histone) %in% levels(aii_merged$histone.y),"color"], 
                   bty="n", cex=0.8, pch=1, )
            # ratio dist
            hist(aii_merged$bitscore.xy, breaks = 40, xlim = c(0,5),
                 col = "gray85", border = "white", main = his,xlab = sprintf("bitscore %s / bitscore secondary", his)) 
            abline(v=1, col="gray", lty=2)
            
            # plot evalue
            plot(aii_merged$evalue.x, aii_merged$evalue.y,
                 col=aii_merged$color.y, log="xy",
                 xlab = sprintf("evalue to %s", his),
                 ylab = "evalue to another histone alignment",
                 main = sprintf("evalue %s",his))
            abline(a=0, b=1, lty=2, col="gray")
            legend("bottomleft", 
                   legend = levels(aii_merged$histone.y), 
                   col = class_col[as.character(class_col$histone) %in% levels(aii_merged$histone.y),"color"], 
                   bty="n", cex=0.8, pch=1, )
            # ratio dist
            hist(log2(aii_merged$evalue.xy), breaks = 40, 
                 col = "gray85", border = "white", main = his,xlab = sprintf("evalue %s / evalue secondary", his)) 
            abline(v=0, col="gray", lty=2)
            
            # ratio of ratios
            plot(aii_merged$bitscore.xy, aii_merged$pident.xy, col=aii_merged$color.y, 
                 main="bitscore and pident ratios",
                 xlab=sprintf("bitscore %s / bitscore secondary",his),
                 ylab = sprintf("pident %s / pident secondary", his)
            )
            abline(a=0, b=1, lty=2, col="gray")
            legend("bottomright", 
                   legend = levels(aii_merged$histone.y), 
                   col = class_col[as.character(class_col$histone) %in% levels(aii_merged$histone.y),"color"], 
                   bty="n", cex=0.8, pch=1, )
            
            # ratio of ratios
            plot(aii_merged$bitscore.xy, aii_merged$evalue.x, col=aii_merged$color.y, 
                 log="y",
                 main="bitscore ratios v evalue",
                 xlab=sprintf("bitscore %s / bitscore secondary",his),
                 ylab = sprintf("evalue %s", his)
            )
            abline(v=1, lty=2, col="gray")
            legend("bottomleft", 
                   legend = levels(aii_merged$histone.y), 
                   col = class_col[as.character(class_col$histone) %in% levels(aii_merged$histone.y),"color"], 
                   bty="n", cex=0.8, pch=1, )
            
            
        } else {
            pie(2)
            pie(3)
            pie(4)
            pie(5)
            pie(6)
            pie(7)
            pie(8)
            pie(9)
        }
    }
    dev.off()
    
    # aii_merged$best_histone = character(length = nrow(aii_merged))
    # aii_merged[aii_merged$bitscore.xy>1 & aii_merged$evalue.x < 1e-30,"best_histone"] = 
    #     as.character(aii_merged[aii_merged$bitscore.xy>1 & aii_merged$evalue.x < 1e-30,"histone.x"])
    
    # classified histones
    # cla = ali[,c("qseqid","histone")]
    # cla = cla[!duplicated(cla[,c("qseqid","histone")]),]
    # 
    # cla = aggregate(ali[,c("qseqid","bitscore")], by=list(ali$qseqid), FUN = sum)
    # 
    # cla = merge(cla, ali, by="qseqid", all.y = F, all.x = T)
    
    # bitscore matrix
    ali_bnr = ali[,c("qseqid","histone","bitscore")]  
    ali_bnr = ali_bnr[!duplicated(ali_bnr[,c("qseqid","histone")]),]
    ali_bwr = reshape(ali_bnr, idvar = "qseqid", direction = "wide", timevar = "histone")
    rownames(ali_bwr) = ali_bwr$qseqid
    ali_bwr = as.matrix(ali_bwr[,-1])
    ali_bwr[is.na(ali_bwr)] = 0 
    
    # pident matrix
    ali_pnr = ali[,c("qseqid","histone","pident")]  
    ali_pnr = ali_pnr[!duplicated(ali_pnr[,c("qseqid","histone")]),]
    ali_pwr = reshape(ali_pnr, idvar = "qseqid", direction = "wide", timevar = "histone")
    rownames(ali_pwr) = ali_pwr$qseqid
    ali_pwr = as.matrix(ali_pwr[,-1])
    ali_pwr[is.na(ali_pwr)] = 0 
    
    cla = data.frame(qseqid = rownames(ali_pwr))
    # cla$length = grep(pattern = "[0-9]+-[0-9]+$",cla$qseqid, value = T)
    
    
    cla$bs_ratio = apply(ali_bwr, 1, function(x) sort(x, decreasing = TRUE)[1] / sort(x, decreasing = TRUE)[2] )
    cla$hist_best = apply(ali_bwr, 1, function(x) names(sort(x, decreasing = TRUE)[1]) )
    cla$hist_best_bs = apply(ali_bwr, 1, function(x) sort(x, decreasing = TRUE)[1] )
    cla$hist_scnd = apply(ali_bwr, 1, function(x) names(sort(x, decreasing = TRUE)[2]) )
    cla$hist_scnd_bs = apply(ali_bwr, 1, function(x) sort(x, decreasing = TRUE)[2] )
    
    cla$classification = "unknown"
    cla[cla$bs_ratio > 1.1,]$classification = cla[cla$bs_ratio > 1.1,]$hist_best
    
    # # both matrices
    # ali_dat = cbind(ali_bwr, ali_pwr)
    # ali_dat_s = scale(ali_dat)
    # 
    # pdf(sprintf("%s.diagnostic_pca.pdf",ali_fi), width = 12, height = 8)
    # layout(mat = matrix(c(1:6),nrow = 2))
    # # PCA
    # message(paste("# PCA"))
    # pic = prcomp(ali_dat_s)
    # pic$variancefraction = pic$sdev^2/sum(pic$sdev^2)
    # 
    # class_colors = rainbow(n=length(levels(as.factor(cla$classification))), v=0.9, end = 0.8)
    # 
    # # Plot PCA: 12
    # for (i in seq(4)) {
    #     for (j in seq(4)) {
    #         if (i < j) {
    #             plot(pic$x[,c(i,j)],col=class_colors[as.factor(cla$classification)],main=sprintf("PCA %i&%i", i,j),
    #                  xlab=sprintf("PCA %i, %.2f%% variance", i,pic$variancefraction[i]*100),
    #                  ylab=sprintf("PCA %i, %.2f%% variance", j,pic$variancefraction[j]*100))
    #             legend("bottomright",legend = levels(as.factor(cla$classification)), col = class_colors[as.factor(levels(as.factor(cla$classification)))],
    #                    pch = 1,cex=0.5)
    #         }
    #         
    #     }
    # }
    # # Plot PCA: eigenvalues
    # plot(pic$variancefraction,type="b",col="red",main="PCA variance")
    # text(pic$variancefraction,labels=signif(pic$variancefrac,digits=3),pos=3,col="black")
    # dev.off()
    # 
    # write.table(cla, file=sprintf("%s.clas.csv", ali_fi), sep = "\t", quote = F, row.names = F)
    
}
