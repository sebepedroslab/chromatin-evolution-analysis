# load libraries
library(ape)
library(scales)
library(phytools)
library(stringr)
library(vioplot)


#### Define input ####

domlist=c("SET","BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","HIF-1","Histone",
          "LinkerHistone","PHD","PTIP","SIR2","SNF2_N","TUDOR",
          "zf-CCHH","zf-CXXC","Hist_deacetyl","Acetyltransf_1")

domlist=c("SET","BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","HIF-1","Histone",
          "LinkerHistone","PHD","PTIP","SIR2","SNF2_N","TUDOR",
          "zf-CCHH","zf-CXXC","Hist_deacetyl")

domlist=c("BIR","Bromodomain","Chromo","CupinJmjC",
          "DOT1","GNAT_acetyltr_2","Hist_deacetyl",
          "Histone","LinkerHistone","Kelch","PHD","PTIP","SAM","SET")

pdf("results_viruses/viral_toolkit_similarity.pdf", width = 10, height = 10)
for (dom in domlist) {
    
    print(dom)
    
    for (taxi in c("euk")) {
        layout(matrix(1:4, nrow = 2))
        dist_ppos = list()
        dist_pide = list()
        for (taxj in c("vir","euk","arc","bac")) {
            if (taxi != taxj) {
                ali = read.table(
                    gzfile(sprintf("results_viruses/pairwise_alignments/diamond.%s.%s-%s.csv.gz", dom, taxj, taxi)),
                    header = F, 
                    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","ppos","qcovhsp"))
                ali = ali[ali$qseqid != ali$sseqid,]
                
                # plot(ali$pident,ali$ppos, col="blue", xlim = c(0,100), ylim = c(0,100))
                # title(sprintf("%s: %s-%s", dom ,taxi, taxj))
                # abline(a=0, b=1, lty=2)
                dist_ppos[[sprintf("%s-%s", taxj, taxi)]] = ali$ppos
                dist_pide[[sprintf("%s-%s", taxj, taxi)]] = ali$pident
            }
        }
        vioplot(dist_ppos, main = sprintf("ppos %s", dom), col = "lightgray", ylim=c(0,100))
        vioplot(dist_pide, main = sprintf("pident %s", dom), col = "lightgray", ylim=c(0,100))
        plot(0,0, xlim = c(0,100), ylim = c(0,100), xlab = "pident", ylab = "ppos")
        points(dist_pide$`arc-euk`, dist_ppos$`arc-euk`, col=alpha("darkturquoise",0.3))
        points(dist_pide$`bac-euk`, dist_ppos$`bac-euk`, col=alpha("blue",0.3))
        points(dist_pide$`vir-euk`, dist_ppos$`vir-euk`, col="deeppink")
        abline(a=0, b=1, lty=2)
        legend("bottomright", col=c("slategray4","blue","darkturquoise","deeppink"), pch=1, legend = c("eukaryotes","bacteria","archaea","viruses"), cex=0.5, bty = "n")
        
    }
    
}
dev.off()
