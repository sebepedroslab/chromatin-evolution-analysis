# libraries
library(data.table)

# input
hido_fn = "histone_domain_coordinates.csv"
clas_fn = "../histone_classification/all.Histone.to_histdb.Spring-Naive.csv"
dole_fn = "domain_lengths.csv"

# parameters
ali_fraction = 0.5

# load domain alignment table
hido = read.table(
  hido_fn, header = F, stringsAsFactors = F, 
  col.names = c("gene","start","end","domains"))
hido$domain_id = paste(hido$gene, paste(hido$start, hido$end, sep="-"), sep="_")
# load classification
clas = read.table(clas_fn, stringsAsFactors = F, header = T, sep = "\t")
clas[clas$classification == "" ,"classification"]  = NA
# load expected domain lengths
dole = read.table(dole_fn, header = F, col.names = c("HMM","HMM_length"))

# find alignment length
hido$length = hido$end - hido$start + 1  # subtract X because on average we've extended domains by this much?

# expected HMM lengths
hido$HMM_length = dole[dole$HMM == "Histone",]$HMM_length # expected HMM lengths: by default, histone
hido[grepl("CENP-T_C", hido$domains),]$HMM_length = dole[dole$HMM == "CENP-T_C",]$HMM_length # if it looks archaeal, archaeal
hido[grepl("CBFD_NFYB_HMF", hido$domains),]$HMM_length = dole[dole$HMM == "CBFD_NFYB_HMF",]$HMM_length # if it looks centromeric, centromeric

# is it a complete domain?
hido$is_complete = hido$length / hido$HMM_length > ali_fraction
hido = hido[hido$is_complete,]

# find candidate double histones
hido_dups_bool = duplicated(hido$gene)
hido_dups_gene = hido[hido_dups_bool, "gene"]
hido_dups = hido[hido$gene %in% hido_dups_gene,]

# add histone type classification
hido_dups = merge(hido_dups, clas, by.x="domain_id", by.y="members", all.x = T, all.y = F)

# reorder domains by coordinate
hido_dups = hido_dups[order(hido_dups$gene, hido_dups$start),]

# add count
hido_dups = data.table(hido_dups)
hido_dups[ , Index := 1:.N , by = c("gene")]

# separate first and second domains
hido_dups_first = hido_dups[hido_dups$Index == 1,]
hido_dups_secnd = hido_dups[hido_dups$Index == 2,]
colnames(hido_dups_first) = paste(colnames(hido_dups_first), "i", sep="_")
colnames(hido_dups_secnd) = paste(colnames(hido_dups_secnd), "j", sep="_")
# rejoin
hido_dups_pairs = cbind(hido_dups_first, hido_dups_secnd)

### archaea double histones ###

# load archaea taxonomy
taxne_fn = "../../data/taxonomy_ranked_wg.tsv.gz"
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=F, stringsAsFactors=F)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")

# load sequence taxonomy (species of origin of each seq)
gentax_fn = "../../data/seq_Archaea.taxa.csv.gz"
gentax = fread(cmd=sprintf("zcat %s", gentax_fn), data.table = F, header = F, stringsAsFactors = F)
colnames(gentax) = c("gene", "taxon")
gentax$gene = paste("arc",gentax$gene, sep="_")
gentax = gentax[gentax$gene %in% hido$gene,]

# get archaeal double histones
hido_dups_pairs = merge(hido_dups_pairs, gentax, all.x = T, all.y = F, by.x="gene_i", by.y = "gene")
hido_dups_pairs = merge(hido_dups_pairs, taxne, all.x = T, all.y = F, by.x="taxon", by.y = "taxon")

# save
write.table(hido_dups_pairs, file="histone_dimers.csv", sep="\t", row.names = F, col.names = T, quote=F)

print(table(hido_dups_pairs$phylum))
