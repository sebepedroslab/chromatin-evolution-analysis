library(data.table)
library(pheatmap)
library(stringr)


# input files
taxne_fn = "../data/taxonomy_ranked_wg.tsv.gz"                 # taxonomic info for non-euks
gen_list_fn = "../data/gene_families_hmm.csv"               # list of gene familes
search_fn = "gene_counts/"
ne_taxon_fo = "../data/"
outp_fn = "results_counts/"
source("../scripts/venn_diagrams.R")


graphics.off()

# heatmap colors
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))
# col_blue = colorRampPalette(interpolate="l",c("gray95", "deepskyblue","dodgerblue3","dodgerblue4"))

# load list hmms
gen_list = read.table(gen_list_fn, header = F, stringsAsFactors = T)
colnames(gen_list) = c("Class","Type","Family","Domains","search_thr","inflation","min_phylo_size")
gen_list_gap_ixs = c(1,1+which(diff(as.numeric(gen_list$Type))!=0)) - 1

#### non-euks HMM counts ####

# create counts tables: for each supergroup (Archaea, Bacteria, Viruses)

# load general taxonomy (species, genus, etc)
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=F, stringsAsFactors=F)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")



#### correlation with histone presence ####


# function: jaccard
jaccard = function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}



noneuk_report_list = list(
    # commented-out taxa have 0 hits
    list(tax="Archaea", tbi="arc", tai="Candidatus Lokiarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Thorarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Odinarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Heimdallarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Crenarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Verstraetearchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Nezhaarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Bathyarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Thaumarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Geothermarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Korarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Diapherotrites", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Micrarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Altiarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Hydrothermarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Aenigmarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Nanoarchaeota", tti="phylum", tfi="genus"),
    list(tax="Archaea", tbi="arc", tai="Candidatus Woesearchaeota", tti="phylum", tfi="genus"),
    list(tax="Bacteria", tbi="bac", tai="Actinobacteria", tti="phylum", tfi="genus"),
    list(tax="Bacteria", tbi="bac", tai="Bacteroidetes", tti="phylum", tfi="genus"),
    list(tax="Bacteria", tbi="bac", tai="Firmicutes", tti="phylum", tfi="genus"),
    list(tax="Bacteria", tbi="bac", tai="Proteobacteria", tti="phylum", tfi="genus"),
    list(tax="Bacteria", tbi="bac", tai="Spirochaetes", tti="phylum", tfi="genus")
)




# list of taxa groups to compare, at TAXON level
noneuk_report_list = list(
    # commented-out taxa have 0 hits
    list(tax="Archaea", tbi="arc", tai="Archaea", tti="superkingdom", tfi="taxon"),
    list(tax="Bacteria", tbi="bac", tai="Actinobacteria", tti="phylum", tfi="taxon"),
    list(tax="Bacteria", tbi="bac", tai="Bacteroidetes", tti="phylum", tfi="taxon"),
    list(tax="Bacteria", tbi="bac", tai="Firmicutes", tti="phylum", tfi="taxon"),
    list(tax="Bacteria", tbi="bac", tai="Proteobacteria", tti="phylum", tfi="taxon"),
    list(tax="Bacteria", tbi="bac", tai="Spirochaetes", tti="phylum", tfi="taxon")
)


# list of domains
dom_test_list = c("Acetyltransf_1","GNAT_acetyltr_2","Hist_deacetyl","SIR2","DOT1","SET","CupinJmjC","ING","MBT","PWWP","SNF2_N","ASF1_hist_chap")
dom_test_mat = matrix(ncol = length(dom_test_list), nrow = 0)
colnames(dom_test_mat) = dom_test_list
rowname_vect = c()


# loop through lists

pdf("results_counts/correlation_with_histones.noneuk-taxon-venn.pdf", width = 2.5, height = 2.5)
for (rep in noneuk_report_list) {
    
    # define subset to plot
    tax=rep$tax # which file to load? (archaea, virus, bacteria)
    tbi=rep$tbi # brief name for the loaded dataset (arc, vir, bac)
    tai=rep$tai # subset to this taxonomic range
    tti=rep$tti # the subset taxonomic range belongs to this category (kingdom...)
    tfi=rep$tfi # focus on reporting at this lower taxonomic level
    
    print(sprintf("gene counts %s (%s) at %s level", tai, tti, tfi))
    
    # load sequence taxonomy (species of origin of each seq)
    gentax_fn = sprintf("%s/seq_%s.taxa.csv.gz", ne_taxon_fo, tax)
    gentax = fread(cmd=sprintf("zcat %s", gentax_fn), data.table = F, header = F)
    colnames(gentax) = c("gene", "taxon")
    
    # restrict taxonomy to group of interest
    taxni = taxne[!is.na(taxne[,tti]),]
    taxni = taxni[taxni[,tti] == tai,]
    # taxni = droplevels(taxni)
    taxni_counts = table(taxni[,tfi])
    taxni_counts = c(sum(taxni_counts),taxni_counts)
    names(taxni_counts)[1] = "Sum"
    print(sprintf("gene counts across %i %s %s(s)", length(unique(taxni[,tfi])), tai, tfi))
    
    # load HMM hits
    gen = fread(input = sprintf("gene_counts/%s_genecounts.csv", tbi), data.table = F, header = F, col.names = c("gene","gene_family"))
    
    # add species info
    gen = merge(gen, gentax, by.x = "gene", by.y = "gene", all.x = T)
    # one row per species
    gen_u = within(gen, rm("gene"))
    gen_u = unique(gen_u)
    
    # add taxonomic info
    gen_wtax = merge(gen_u, taxni, by.x = "taxon", by.y = "taxon", all.x = T, all.y = F)
    
    # taxonomic groupings and gene families as factors
    levels_taxon_factor = levels(as.factor(taxni[,tfi]))
    gen_wtax[,"taxon_factor"] = factor(gen_wtax[,tfi], levels = levels_taxon_factor)
    gen_wtax[,"genefam_factor"] = factor(gen_wtax$gene_family, levels = as.vector(gen_list$Family))
    
    # crosstabulation
    gen_crosstab = xtabs(formula = ~ genefam_factor + taxon_factor, data = gen_wtax, drop.unused.levels	= F)
    
    # binarise (presence per taxon)
    gen_crosstab [ gen_crosstab > 1 ] = 1
    
    
    # jaccard index between histone presence and each domain
    hjac_tai = unlist(lapply(dom_test_list, function(dom) { 
        tax_with_hist = colnames(gen_crosstab)[ gen_crosstab["Histone",]>0 ]
        tax_with_domi = colnames(gen_crosstab)[ gen_crosstab[dom,]>0 ]
        jaccard(tax_with_hist, tax_with_domi)
    }))
    names(hjac_tai) = dom_test_list
    
    # store
    rowname_vect = c(rowname_vect, paste(tai, tfi))
    dom_test_mat = rbind(dom_test_mat, hjac_tai)
    
    # venn diagrams
    for (dom in dom_test_list) {
        tax_with_hist = colnames(gen_crosstab)[ gen_crosstab["Histone",]>0 ]
        tax_with_domi = colnames(gen_crosstab)[ gen_crosstab[dom,]>0 ]
        venn.two(tax_with_hist, tax_with_domi, catname1 = "Histone", catname2 = dom, 
                 main = sprintf("%s | %s %s\nJaccard = %.2f", dom, tai, tfi, hjac_tai[dom]))
    }
    
}
dev.off()

# plot at taxon level
pdf("results_counts/correlation_with_histones.noneuk-taxon.pdf", width = 5, height = 3)
rownames(dom_test_mat) = rowname_vect
col_tax = rainbow(n=nrow(dom_test_mat), start = 0.1, end = 0.85, v = 0.9)
barplot(dom_test_mat, las=2, beside = T, col = col_tax, cex.names=0.6, cex.axis = 0.6, ylab = "Jaccard index", cex.lab=0.6)
title(main="Jaccard index with Histone pres/abs profile", cex.main=0.6)
legend("topright", legend = rownames(dom_test_mat), fill = col_tax, cex=0.4)

write.table(dom_test_mat, file="results_counts/correlation_with_histones.noneuk-taxon.csv", row.names = T, col.names = T, quote = F, sep = "\t")
dev.off()

