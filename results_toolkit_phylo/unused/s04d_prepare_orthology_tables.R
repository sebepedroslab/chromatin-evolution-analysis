#### Input ####
# load libraries
library(stringr)

dat_fn  = "orthogroups_euk.csv"
phyl_fn = "data/species_tree.newick"
clas_fn = "../data/gene_families_hmm.csv"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"


#### Define input ####

# load data
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

# read species tree
phyl = ape::read.tree(file = phyl_fn)
phyl$edge.length = rep(1, nrow(phyl$edge))
# dataframe of edges
phyl_edge           = as.data.frame(phyl$edge)
colnames(phyl_edge) = c("edge_start","edge_end")
phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
# dataframe of nodes
phyl_nods = data.frame(taxa = c(phyl$tip.label, phyl$node.label))
phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
phyl_nods$is_tip   = phyl_nods$edge_end <= length(phyl$tip.label)

# merge them
phyl_edge = merge(phyl_edge, phyl_nods, all.x = T, all.y = T, by.x = "edge_end", by.y = "edge_end")



#### Table with all genes ####

# add gene family
out = dat
out[,"Gene family"] = stringr::str_split(out$orthogroup, "\\.", simplify = T)[,1]
out[,"sps"] = stringr::str_split(out$gene, "_", simplify = T)[,1]

out = merge(out, cla, by.x = "Gene family", by.y = "gene_fam")
out = merge(out, tax, by.x = "sps", by.y = "Species")
out = out[, c("gene","Species.name","gene_type","Gene family","orthogroup","orthologous_to")]
colnames(out) = c("Gene","Species","Category","Gene family","Orthogroup","Direct orthologs in human")

out$Species = factor(out$Species, levels = tax$Species.name)

out = out[order(out$Orthogroup, out$Species),]


write.table(out, "data/orthogroups_euk_annotated.csv", sep="\t", quote = F, row.names = F, col.names = T)


#### Table with reference species ####

ref_sps = c("Hsap","Dmel","Scer","Atha")
dar = dat
dar$sps = stringr::str_split(dar$gene, "_", simplify = T)[,1]
dar = dar [ dar$sps %in% ref_sps, c("gene","orthogroup","sps")]

# get gene symbols
greg = data.frame()
for (spi in ref_sps) { 
    
    gref = read.table(sprintf("data/gene_names_%s.csv", spi), header = F, col.names = c("gene","Gene symbol"))
    greg = rbind(greg, gref)
    
}

dar = merge(dar, greg,by.x = "gene", by.y = "gene", all.x = T, all.y = F)
dar$sps = factor(dar$sps, levels = ref_sps)

# concatenate genes from the same OG
dagg = aggregate(dar$gene, by = list(Orthogroup = dar$orthogroup, sps = dar$sps), FUN = function(x) paste(x[!is.na(x)], collapse = ", "))
colnames(dagg)[3] = "Genes"

# concatenate gene symbols from the same OG
dags = aggregate(dar$Gene.symbol, by = list(Orthogroup = dar$orthogroup, sps = dar$sps), FUN = function(x) paste(x[!is.na(x)], collapse = ", "))
colnames(dags)[3] = "Gene symbols"

# paste
daga = cbind(dagg, dags$`Gene symbols`)
colnames(daga)[4] = "Gene symbols"

# widen to species
daw = reshape(daga, idvar = "Orthogroup", timevar = "sps", direction = "wide")
daw = daw[order(daw$Orthogroup),]

for (col in 2:ncol(daw)) {
    daw [is.na(daw[,col]) , col]  = ""
}

colnames(daw) = gsub("\\."," ",colnames(daw))
colnames(daw) = gsub("Atha","Arabidopsis thaliana",colnames(daw))
colnames(daw) = gsub("Scer","Saccharomyces cerevisiae",colnames(daw))
colnames(daw) = gsub("Dmel","Drosophila melanogaster",colnames(daw))
colnames(daw) = gsub("Hsap","Homo sapiens",colnames(daw))

write.table(daw, "data/orthogroups_euk_ref_gene_symbols.csv", sep="\t", quote = F, row.names = F, col.names = T)
