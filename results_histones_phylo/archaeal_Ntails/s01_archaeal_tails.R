# libraries
library(stringr)
library(seqinr)
library(data.table)

# input files
cla_fn = "../histone_classification/all.Histone.to_histdb.Spring-Naive.csv"
pfa_fn = "arc.Histone.seqs.pfamscan.csv"
fas_fn = "arc.Histone.seqs.fasta"
vtax_fn = "../../data/seq_Archaea.taxa.csv.gz"    # taxon for each viral seq
taxne_fn = "../../data/taxonomy_ranked_wg.tsv.gz" # taxonomy for non-euks


# read data
# taxonomy of viral sequences
vtax = fread(cmd=sprintf("zcat %s", vtax_fn), data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
colnames(vtax) = c("gene", "taxon")
# taxonomy of species
taxne = fread(cmd=sprintf("zcat %s", taxne_fn), sep = "\t", na.strings = "", quote = "", data.table=FALSE, stringsAsFactors=FALSE)
colnames(taxne) = c("taxnid","taxon","species","genus","family","order","class","phylum","kingdom", "superkingdom")
taxne$phylum = gsub("Candidatus ", "",taxne$phylum)


# read fasta
fas = seqinr::read.fasta(fas_fn)
names(fas) = stringr::str_split(names(fas), " ", simplify = TRUE)[,1]

# read tables: PFAM searches
pfa = read.table(
  pfa_fn, stringsAsFactors = FALSE,
  col.names = c("gene","ali_start","ali_end","env_start","env_end","dom_acc","dom","dom_type","hmm_start","hmm_end","hmm_len","bitscore","eval","significant","dom_clan"))
pfa$ali_len = pfa$ali_end - pfa$ali_start
pfa$env_len = pfa$env_end - pfa$env_start
# pfa$seq_id  = paste("arc_",pfa$gene,"_",pfa$ali_start,"-",pfa$ali_end, sep = "")

# add taxonomy
pfa = merge(pfa, vtax, by.x = "gene", by.y = "gene", all.x = TRUE, all.y=FALSE)
pfa = merge(pfa, taxne, by.x = "taxon", by.y = "taxon", all.x = TRUE, all.y = FALSE)

# classification
cla = read.table(cla_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cla$gene = gsub("^arc_","", cla$members)
cla$gene = gsub("_\\d+-\\d+$","", cla$gene)

# keep only genes with domain of interest, a minimum fraction of the original HMM model, and belonging to the archaeal component
pfi = pfa [ pfa$dom == "CBFD_NFYB_HMF", ]
# pfi = pfi [ pfi$env_len > pfi$hmm_len * 0.7, ]
pfi = pfi [ pfi$hmm_start < 10 & pfi$hmm_end > 55 , ]
pfi = merge(pfi, cla, by.x = "gene", by.y = "gene", all.x = TRUE, all.y = FALSE)
pfi = pfi [ pfi$component == 0  & !is.na(pfi$component), ]


# get sequences that have a complete domain and at least 10 aa BEFORE this domain
pfo = pfi [ pfi$env_start > 10, ]

# get sequences in correct order
fai = fas[pfo$gene]

# get core and tail sequence
pfo$seq_tail = unlist(lapply(1:length(fai), function(x) paste(fai[[x]] [ 1:pfo$env_start[x] - 1 ] , collapse = "")))
pfo$seq_core = unlist(lapply(1:length(fai), function(x) paste(fai[[x]] [ pfo$env_start[x]:pfo$env_end[x] ] , collapse = "")))
pfo$seq_tail_len = unlist(lapply(1:length(fai), function(x) length(fai[[x]] [ 1:pfo$env_start[x] - 1 ])))
pfo$seq_core_len = unlist(lapply(1:length(fai), function(x) length(fai[[x]] [ pfo$env_start[x]:pfo$env_end[x] ])))

# get upper case
pfo$seq_core = toupper(pfo$seq_core)
pfo$seq_tail = toupper(pfo$seq_tail)


# get K freq
pfo$seq_tail_Kfreq = unlist(lapply(1:length(fai), function(x) sum(unlist(strsplit(pfo$seq_tail[x],"")) == "K") / pfo$seq_tail_len[x]  ))
# order
pfo = pfo [order(pfo$seq_tail_Kfreq, decreasing = TRUE) ,]

# remove duplicates
pfo = pfo [ !duplicated(pfo$gene), ]

# save plot 
pdf("summary_archaeal_tails_Kfreq.pdf", height = 5, width = 3.5)

pfo_phylum = pfo$phylum
pfo_phylum [ is.na(pfo_phylum) ] = "none"
cols = rainbow(n=nlevels(as.factor(pfo_phylum)), v = 0.8, start = 0.1, end = 0.9)
names(cols) = levels(as.factor(pfo_phylum))

# plot
plot(pfo$seq_tail_Kfreq, col=cols[pfo_phylum], ylab = "K frequency", pch=1,las=1, cex=sqrt(pfo$seq_tail_len / 15))
title(main=sprintf("K frequency n=%i", nrow(pfo)))
legend("topright",legend=names(cols), col=cols, pch=1, bty = "n", cex=0.5)

# only from selected taxon
tax_list = c("Methanobrevibacter cuticularis", "Nitrososphaera viennensis", "Methanosarcina spelaei", "Methanospirillum stamsii")
tax_bool = pfo$taxon %in% tax_list
pfo_i = pfo[tax_bool,]
plot(x=which(tax_bool),y=pfo$seq_tail_Kfreq[tax_bool], col=cols[pfo_phylum][tax_bool], ylab = "K frequency", pch=1,las=1, cex=sqrt(pfo$seq_tail_len / 15)[tax_bool],
		 ylim = c(min(pfo$seq_tail_Kfreq),max(pfo$seq_tail_Kfreq)), xlim = c(1,84))
title(main=sprintf("K frequency n=%i", nrow(pfo[tax_bool,])))
legend("topright",legend=names(cols), col=cols, pch=1, bty = "n", cex=0.5)


dev.off()

# save table
write.table(pfo, file = "summary_archaeal_tails.csv", quote = FALSE, sep = "\t", row.names = FALSE)

# save tails in FASTA
seqinr::write.fasta(as.list(pfo$seq_tail), as.string = TRUE, names = paste(pfo$gene, pfo$taxon), file.out = "tails_archaea.fasta")

