# input
cla_fn = "../data/gene_families_hmm.csv"
tax_fn = "../data/euk_taxonomy_annotated_2020-08-11.csv"
dat_fn = "orthogroups_euk.csv"
# ogs_fn = "orthogroups_euk.ancestral.counts.csv"
arq_fn = "architectures.csv"
# ancestral presence probabilities
# ogp_fn = "orthogroups_euk.ancestral.posteriors_pres.csv"
# ogw_fn = "orthogroups_euk.ancestral.wagner_g5_pres.csv"
pho_fn = "results_preeuk/summary_preeuk_families.annotation_per_OG.csv"
ane_fn = "results_evol/leca_reconstruction_evidence.csv"

#### Load ####

# read
cla = read.table(cla_fn, header = FALSE, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"), stringsAsFactors = TRUE)
tax = read.table(tax_fn, sep = "\t", header = TRUE, stringsAsFactors = TRUE)
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE)
arq = read.table(arq_fn, sep = "\t", stringsAsFactors = FALSE, header = FALSE, col.names = c("gene","architecture"))

# factor
cla$gene_type = factor(
  cla$gene_type, 
  levels = c("Acetylase","Deacetylase","Methylase","Demethylase","Readers","Remodeller","Chaperones","Histones","KMT1","KMT2","KMT4","PC1","PC2","common") )

# load ancestral reconstruction
# ogp = read.table(ogp_fn)
# ogw = read.table(ogw_fn, comment.char =  "", header = TRUE)
# rownames(ogw) = gsub(":.*","", ogw$name)
# ogw = ogw[,4:ncol(ogw)]
ane = read.table(ane_fn, header = TRUE, sep = "\t")
pho = read.table(pho_fn, header = TRUE)


#### Init classification ####

# report at OG level
pho$gene_fam = gsub("\\.HG.*","", pho$orthogroup)
out = merge(pho, cla[,c("gene_fam","gene_class","gene_type")], by.x = "gene_fam", by.y = "gene_fam", all.x = TRUE, all.y = FALSE)
out$orthogroup_short = gsub(":.*","",out$orthogroup)
out$closest_group [ out$closest_group == "eukoth" ] = "Other eukaryotic"
out$closest_group [ out$closest_group == "bac" ] = "Bacteria"
out$closest_group [ out$closest_group == "arc" ] = "Archaea"
out$closest_group [ out$closest_group == "vir" ] = "Virus"

# add ancestral LECA support
out = merge(out, ane[,c("og_id", "ancestral_support_string")], by.x = "orthogroup_short", by.y = "og_id", all.x = TRUE, ally = FALSE)


#### Presence count summary ####

# count homologs, num sps, num lineages... per og
dat$orthogroup_short = gsub(":.*","",dat$orthogroup)
dat_a = data.frame(table(dat$orthogroup_short))
colnames(dat_a) = c("orthogroup_short", "num_genes")

# num of species & macrotaxa per OG
dat$species = gsub("_.*","", dat$gene)
dat = merge(dat, tax[,c("Species","Macrogroup")], by.x = "species", by.y = "Species")
dat_a$num_species = aggregate(species ~ orthogroup_short, data = dat, function(i) length(unique(i)) )[,2]
dat_a$num_lineages = aggregate(Macrogroup ~ orthogroup_short, data = dat, function(i) length(unique(i)) )[,2]
dat_a$lineages_string = aggregate(Macrogroup ~ orthogroup_short, data = dat, function(i) paste(unique(i), collapse = ",") )[,2]

# merge with presence count summary
out = merge(out, dat_a, by = "orthogroup_short", all.x = TRUE, all.y = FALSE)


#### Model species ####

# other model species
mod_a = read.table("data/gene_names_Atha.csv", col.names = c("gene","Atha_symbol"))
mod_y = read.table("data/gene_names_Scer.csv", col.names = c("gene","Scer_symbol"))
mod_d = read.table("data/gene_names_Dmel.csv", col.names = c("gene","Dmel_symbol"))
mod_h = read.table("data/gene_names_Hsap.csv", col.names = c("gene","Hsap_symbol"))

mod_a = merge(mod_a, dat, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
mod_y = merge(mod_y, dat, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
mod_d = merge(mod_d, dat, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
mod_h = merge(mod_h, dat, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
mod_a_a = aggregate(Atha_symbol ~ orthogroup_short, data = mod_a, function(i) paste(unique(i), collapse = ",") )
mod_y_a = aggregate(Scer_symbol ~ orthogroup_short, data = mod_y, function(i) paste(unique(i), collapse = ",") )
mod_d_a = aggregate(Dmel_symbol ~ orthogroup_short, data = mod_d, function(i) paste(unique(i), collapse = ",") )
mod_h_a = aggregate(Hsap_symbol ~ orthogroup_short, data = mod_h, function(i) paste(unique(i), collapse = ",") )

# add to output
out = merge(out, mod_h_a, all.x = TRUE, all.y = FALSE, by.x = "orthogroup_short", by.y = "orthogroup_short")
out = merge(out, mod_d_a, all.x = TRUE, all.y = FALSE, by.x = "orthogroup_short", by.y = "orthogroup_short")
out = merge(out, mod_y_a, all.x = TRUE, all.y = FALSE, by.x = "orthogroup_short", by.y = "orthogroup_short")
out = merge(out, mod_a_a, all.x = TRUE, all.y = FALSE, by.x = "orthogroup_short", by.y = "orthogroup_short")



#### Consensus architecture ####

# load
gcla_fn = "gene_counts/euk_genecounts.csv"
arq = read.table(arq_fn, header = F, sep="\t", col.names = c("gene","architecture"))
gcla = read.table(gcla_fn, header = F, sep = "\t", col.names = c("gene","gene_fam"))
arq_cla = merge(arq, gcla, by.x = "gene", by.y="gene", all.x = F, all.y = T)
arq_cla = merge(arq_cla, cla, by.x = "gene_fam", by.y = "gene_fam", all.x = T, all.y = F)

# super-simplify architectures (collapse consecutive repeats)
simple_arqs_list = stringr::str_split(arq_cla$architecture, pattern = " ")
# simple_arqs_list = sapply(simple_arqs_list, function(x) gsub(pattern = "_\\d+$", "", x))
simple_arqs_vect = sapply(simple_arqs_list, FUN = function(x) paste(unique(x), collapse = ",") )
arq_cla$simarq = simple_arqs_vect
# add ogs
arq_cla = merge(arq_cla, dat, by.x = "gene", by.y = "gene", all.x = TRUE)

arq_a = aggregate(simarq ~ orthogroup_short, data = arq_cla, function(i) {
  tf = sort(table(i), decreasing = TRUE) / sum(table(i))
  if (tf[1] >= 0.25) {
    st = names(tf)[1]  
  } else {
    st = "No consensus"
  }
  return(st)
})

# add to output table
out = merge(out, arq_a, by = "orthogroup_short", all.x = TRUE, all.y = FALSE)

#### Reformat ####

# reformat output and write file
orthogroup_name_t = stringr::str_split(out$orthogroup, pattern = ":", simplify = TRUE)[,2:3]
orthogroup_name_t = apply(orthogroup_name_t, 1, function(i) paste(i, collapse = ":"))
orthogroup_name_t = gsub(":$","",orthogroup_name_t)
out$orthogroup_name = orthogroup_name_t

# get columns
ouc = out [ c(
  "orthogroup_short",
  "orthogroup_name",
  "gene_fam",
  "gene_type",
  "num_genes",
  "num_species",
  "num_lineages",
  "LECA_probability",
  "ancestral_support_string",
  "closest_group",
  "closest_score",
  "simarq",
  "Hsap_symbol",
  "Dmel_symbol",
  "Scer_symbol",
  "Atha_symbol")]

# rename
colnames(ouc) = c(
  "Orthogroup",
  "Orthogroup annotation",
  "Gene family",
  "Gene type",
  "Num genes",
  "Num species",
  "Num lineages",
  "LECA presence probability",
  "LECA presence support",
  "Closest homologs in...",
  "Closest homologs score",
  "Consensus architecture",
  "Hsap symbol",
  "Dmel symbol",
  "Scer symbol",
  "Atha symbol")

# reorder
ouc = ouc[ order(ouc$`Gene type`,ouc$`Gene family`, -ouc$`Num lineages`, -ouc$`LECA presence probability`),  ]

# reformat floating strings
ouc[,"LECA presence probability"] = sprintf("%.3f", ouc[,"LECA presence probability"])
ouc[,"Closest homologs score"]    = sprintf("%.3f", ouc[,"Closest homologs score"])


# write out
write.table(ouc, file = "orthogroups_euk.summary.csv", sep = "\t", row.names = FALSE, quote = FALSE)





#### annotate per-gene table ####

dat$gene_fam = gsub("\\..*","",dat$orthogroup)
dog = merge(dat, cla, by.x = "gene_fam", by.y = "gene_fam", all.x = TRUE, all.y = FALSE)
dog = merge(dog, tax, by.x = "species", by.y = "Species", all.x = TRUE, all.y = FALSE)

# og name
orthogroup_name_t = stringr::str_split(dog$orthogroup, pattern = ":", simplify = TRUE)[,2:3]
orthogroup_name_t = apply(orthogroup_name_t, 1, function(i) paste(i, collapse = ":"))
orthogroup_name_t = gsub(":$","",orthogroup_name_t)
dog$orthogroup_name = orthogroup_name_t

# reorder
dog = dog [ order( match(dog$orthogroup_short, ouc$Orthogroup), match(dog$species, tax$Species) ) , ]
dog = dog [ , c("gene","Species.name","Group","orthogroup_short","orthogroup_name","gene_fam","gene_type") ]

# write out
write.table(dog, file = "orthogroups_euk.table.csv", sep = "\t", row.names = FALSE, quote = FALSE)

