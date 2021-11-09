# libraries
library("stringr")
library("ape")
library("phangorn")
library("data.table")
library("dbscan")
library("adephylo")
library("vioplot")
library("readxl")

# input
tre_li = list(
    list(file = "alignments/euk.Histone.to_histdb/CC2_H3.iqt.treefile", name="H3", color="blue", eps=1.5, dro=1),
    list(file = "alignments/euk.Histone.to_histdb/CC1_H4.iqt.treefile", name="H4", color="orange", eps=1, dro=1),
    list(file = "alignments/euk.Histone.to_histdb/CC4_H2A.iqt.treefile", name="H2A", color="springgreen3", eps=2, dro=1.5),
    list(file = "alignments/euk.Histone.to_histdb/CC0_H2B.iqt.treefile", name="H2B", color="purple", eps=5, dro=1.5)
)
aln_fn = "euk.Histone.to_histdb.csv.gz"
tax_fn = "../../data/euk_taxonomy_annotated_2020-08-11.csv"

# histones that appear in the proteomics datasets
pro_fn = "all.Histone.to_proteomics.csv.gz"
mod_fn = "../../results_PTMs/consensus_modifications_perseq.xlsx" # list of canonical histones used in proteomics analyse
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


# load ali
aln = read.table(
    gzfile(aln_fn), sep = "\t", 
    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
    stringsAsFactors = FALSE)

# load taxonomy (species order)
tax = read.table(tax_fn, sep="\t", header = TRUE, stringsAsFactors = FALSE)
sps_list = tax$Species
tax$Species_factor = factor(tax$Species, levels = tax$Species)


# define seed canonical sequences as...
pid=100

# prepare outputs
dcn_li = list()
dnn_li = list()
dcc_li = list()
dcns_li = list()

pdf("../trees_histone_distance.pdf", height = 9, width = 16)
layout(matrix(1:8, ncol = 4))
for (tri in tre_li) {
    
    # data
    tre_fn  = tri$file
    his = tri$name
    hco = tri$color
    eps = tri$eps
    dro = tri$dro
    
    # load tree
    message(sprintf("histone %s | load", his))
    tre = ape::read.tree(tre_fn)
    tre = phangorn::midpoint(tre)    
    tre = ape::ladderize(tre) 
    # ix_tip = 1:length(tre$tip.label)
    # ix_tip_edges = which(tre$edge[,2] %in% ix_tip)
    # ix_tip_nodes = tre$edge[,2][ix_tip_edges]
    
    # identify a group of canonical histones (large group of highly similar histones)
    ali = aln[aln$qseqid %in% tre$tip.label,]
    ali_f = ali[ali$pident >= pid & !duplicated(ali$qseqid) & grepl(sprintf("%s$",his), ali$sseqid),]
    canonical_tips_seed = unique(ali_f$qseqid)

    # internode distances    
    # tre_phydis = ape::cophenetic.phylo(tre)
    message(sprintf("histone %s | patristic", his))
    tre_phydis = adephylo::distTips(x=tre, method="patristic")
    tre_phydis = as.matrix(tre_phydis)
    # tre_noddis = ape::node.depth.edgelength(tre)
    # tre_tipdis = tre_noddis[ix_tip_edges]

    # distance to canonical
    mean_cancan = median(tre_phydis[canonical_tips_seed,canonical_tips_seed])
    tre_phydis_tocan = tre_phydis[canonical_tips_seed,]
    closest_can = apply(tre_phydis_tocan, 2, min)
    tips_canonical = names(closest_can[closest_can < mean_cancan])
    
    # which tips are canonical based on distance to a canonical sequence?
    ix_canonical = which(tre$tip.label %in% tips_canonical)
    ix_noncanonical = which(!(tre$tip.label %in% tips_canonical))
    tips_noncanonical = tre$tip.label[ix_noncanonical]
    
    # inter-tip distances
    message(sprintf("histone %s | plot", his))
    hist(tre_phydis, col=hco, breaks = 60, border = "white", xlim=c(0,6), xlab="Distance", main="pairwise distance", sub=sprintf("Median canonical-canonical distance = %.3f", mean_cancan), ylim = c(0,1.4e5), las = 1)
    abline(v=mean_cancan,lty=2, col="black")
    m=0
    for (n in 1:length(hip_list)) {
        hip = hip_list[n]
        hic = hip_cols[n]
        aii_hip = tre_phydis [ ix_canonical , colnames(tre_phydis) %in% his_in_pro[[hip]]  ]
        if (length(aii_hip) > 0 ) {
            # points(x=aii_hip, y = rnorm(length(aii_hip), mean = n * 100, sd = 10), col = alpha(hic,0.6))
            vioplot(aii_hip, at = 1.2e5 + m * 9000, horizontal = TRUE, col = hic, ylim = c(0,100), add = TRUE, wex = 8000)
            m = m + 1
        }
    }
    legend("topright",hip_list, col = hip_cols, pch = 1, bty = "n")
    
    # # cluster with dbscan (assume canonical are large clusters)
    # tre_phyclu = dbscan::dbscan(tre_phydis, eps=eps, minPts = 40)
    # tre_phyclu_df = data.frame(domain = rownames(tre_phydis), cluster = tre_phyclu$cluster)
    # tre_phyclu_top = tre_phyclu_df[tre_phyclu_df$cluster != 0 ,]
    # canonical_tips = as.character(tre_phyclu_top$domain)
    # # which tips are canonical based on dbscan?
    # ix_canonical = which(rownames(tre_phydis) %in% canonical_tips)
    # ix_noncanonical = which(!(colnames(tre_phydis) %in% canonical_tips))
    # 
    # # which tips are canonical based on distance to root?
    # ix_canonical = which(rownames(tre_phydis) %in% canonical_tips)
    # ix_noncanonical = which(!(colnames(tre_phydis) %in% canonical_tips))
    
    # store distance to canonical-noncanonical
    dcn_li[[his]] = apply(tre_phydis [ix_canonical,ix_noncanonical], 2, mean)
    # store distance between noncanonical
    dnn_li[[his]] = apply(tre_phydis [ix_noncanonical,ix_noncanonical], 2, mean)
    # store distance between canonical
    dcc_li[[his]] = apply(tre_phydis [ix_canonical,ix_canonical], 2, mean)
    
    # find average distance between canonical and non canonical histones, for each species
    tip_species = stringr::str_split(rownames(tre_phydis), "_", simplify = TRUE)[,1]
    sps_list_dis = lapply(sps_list, function(spi) { 
        unique(tre_phydis [
            tip_species == spi & rownames(tre_phydis) %in% tips_noncanonical,
            tip_species == spi & rownames(tre_phydis) %in% tips_canonical
        ])
    })
    names(sps_list_dis) = sps_list

    # store distance to canonical-noncanonical
    dcns_li[[his]] = unlist(lapply(sps_list_dis, mean))
    
    # plot tree
    ape::plot.phylo(tre, edge.color = hco, underscore = TRUE, cex = 0.2, font=1, col="grey",type = "p", show.tip.label =  FALSE,main=his)
    ape::add.scale.bar(x=0, y=-10, lcol="darkgray", cex=0.5,length = 0.1)
    ape::tiplabels(pch = 19, col = "red", height = 4, cex = 0.5, tip = ix_canonical)
    
    for (m in 1:length(hip_list)) {
        hip = hip_list[m]
        hic = hip_cols[m]
        ix_hip = which(tre$tip.label %in% his_in_pro[[hip]])
        if (length(ix_hip) > 0) {
            ape::tiplabels(tip = ix_hip, pch=19, col=hic, cex=0.5)
        }
    }

}

# canonical-noncanonical
boxplot(dcn_li, col=c("blue","orange","springgreen3","purple"), outline = TRUE, ylab = "Distance", main="distance canonical-noncanonical")
plot(NA,NA,xlim = c(0,5), ylim = c(0,1), xlab = "Distance", ylab="Fraction")
lines(ecdf(dcn_li$H3), col="blue",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcn_li$H4), col="orange",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcn_li$H2A), col="springgreen3",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcn_li$H2B), col="purple", verticals=TRUE, do.points=FALSE)

# canonical-noncanonical
boxplot(dcns_li, col=c("blue","orange","springgreen3","purple"), outline = TRUE, ylab = "Distance", main="distance canonical-noncanonical, species-wise")
plot(NA,NA,xlim = c(0,5), ylim = c(0,1), xlab = "Distance", ylab="Fraction")
lines(ecdf(dcns_li$H3), col="blue",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcns_li$H4), col="orange",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcns_li$H2A), col="springgreen3",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcns_li$H2B), col="purple", verticals=TRUE, do.points=FALSE)

# noncanonical
boxplot(dnn_li, col=c("blue","orange","springgreen3","purple"), outline = TRUE, ylab = "Distance", main="distance between noncanonical")
plot(NA,NA,xlim = c(0,5), ylim = c(0,1), xlab = "Distance", ylab="Fraction")
lines(ecdf(dnn_li$H3), col="blue",verticals=TRUE, do.points=FALSE)
lines(ecdf(dnn_li$H4), col="orange",verticals=TRUE, do.points=FALSE)
lines(ecdf(dnn_li$H2A), col="springgreen3",verticals=TRUE, do.points=FALSE)
lines(ecdf(dnn_li$H2B), col="purple", verticals=TRUE, do.points=FALSE)

# canonical
boxplot(dcc_li, col=c("blue","orange","springgreen3","purple"), outline = TRUE, ylab = "Distance", main="distance between canonical")
plot(NA,NA,xlim = c(0,2.5), ylim = c(0,1), xlab = "Distance", ylab="Fraction")
lines(ecdf(dcc_li$H3), col="blue",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcc_li$H4), col="orange",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcc_li$H2A), col="springgreen3",verticals=TRUE, do.points=FALSE)
lines(ecdf(dcc_li$H2B), col="purple", verticals=TRUE, do.points=FALSE)



dev.off()

