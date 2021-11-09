# load libraries
library(stringr)
library(GenomicRanges)
library(ape)
library(pheatmap)
library(alluvial)
library(scales)

graphics.off()

#### Define input ####

clas_fn = "../data/gene_families_hmm.csv"
tecl_fn = "../data/Dfam.TEclasses.txt"
arqs_fn = "results_TEfusions/architectures_with_TEs.csv"
tedo_fn = "../data/transposable_element_domains.csv"
alis_fn = "results_TEfusions/all_hits.Dfam.tblastn.csv"
vali_fn = "../data/validation_species.csv"
gecl_fn = "gene_counts/euk_genecounts.csv"
phyl_fn = "../results_toolkit_phylo/data/species_tree.newick"
orth_fn = "../results_toolkit_phylo/orthogroups_euk.csv"
orte_fn = "gene_trees_TE/"
arqs_fo = "results_domains/"


# genes for which data is taken from transcriptome assemblies will be assumed to be supported by full evidence: 
# - expressed, obviously
# - full-length expression and contiguous assembly (even if NNN, bridge should be supported by paired end reads)
sps_transcriptome_only = c("Jaklib","Secula","Andgod","Spimul","Malcal","Maljak","Plamic","Diprot","Gonavo","Masbal","Rigram")

#### Load metadata ####

# load tables
cla = read.table(clas_fn, header = F, sep = "\t", 
                 col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"),
                 stringsAsFactors = F)
# temporarily amend core domain for polycomb RING proteins (include zf-C3HC4_2 in addition to RAWUL)
cla[cla$gene_fam == "PCRing","domains"] = "zf-C3HC4_2,RAWUL"
cla_core_domains = stringr::str_split(paste(cla$domains, collapse = ","), ",")[[1]]

# load various per-gene tables
gec = read.table(gecl_fn, sep = "\t", stringsAsFactors = F, header = F, col.names = c("gene","class"))
val = read.table(vali_fn, sep = "\t", stringsAsFactors = F, header = T)
ali = read.table(alis_fn, sep = "\t", col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), stringsAsFactors = F)
tec = read.table(tecl_fn, sep = "\t", col.names = c("TE_id","TE_family","TE_class"), stringsAsFactors = F)
ted = read.table(tedo_fn, sep= "\t", header = T, stringsAsFactors = F)
ort = read.table(orth_fn, sep = "\t", header = T, stringsAsFactors = F)


# TE taxonomy
TE_class_tax = stringr::str_split(unique(tec$TE_class), pattern = ";", simplify = T)

# orthology family class
ort$family = stringr::str_split(ort$orthogroup, pattern = "\\.", simplify = T)[,1]

# species to report on
val = val[!is.na(val$SRA_runs),]
val = val[val$Num.genes>0,]
val$num_experiments = stringr::str_count(val$SRA_runs,",")+1
val_sps_list = val$Species
# sps_list = c("Drer","Bralan","Spur","Chabra","Rhidel","Rirr")

# load candidate TE fusion genes
can = read.table(arqs_fn, sep = "\t", col.names = c("gene", "orthogroup", "orthologous_to", "gene_family","architecture"), stringsAsFactors = F)
can = within(can, rm("orthogroup","orthologous_to","gene_family"))
can = unique(can)
can$species = stringr::str_split(can$gene, pattern = "_", simplify = T)[,1]
can_sps_list = unique(can$species)

# load TE orthogroups
orte_list = list.files(path = orte_fn, pattern = "*_groups.csv", full.names = T)
# orte_list = list.files(path = orte_fn, pattern = "euk.RVT.*_groups.csv", full.names = T)
orte = data.frame()
for (orti in orte_list) {
  orte = rbind(orte, read.table(orti, header = T, sep = "\t"))
}
# rename domains to genes
colnames(orte) = c("domain","orthogroup_TE","orthologous_to_TE")
orte$gene = gsub(pattern = "_\\d+-\\d+$", replacement="",orte$domain)
orte$fusedTE = stringr::str_split(orte$orthogroup_TE, pattern = "\\.", simplify = T)[,1]
orte$orthogroup_TE = stringr::str_split(orte$orthogroup_TE, pattern = ":", simplify = T)[,1]
# keep only one annotation per gene
# orte = orte[!duplicated(paste(orte$gene, orte$fusedTE)),]
orte_sum = aggregate(orte$orthogroup_TE, by=list(orte$gene), function(x) paste(unique(x), collapse=","))
colnames(orte_sum) = c("gene","orthogroup_TE")

# order species (from phylogeny)
phyl = read.tree(file = phyl_fn)
sps_order = phyl$tip.label
can_sps_list = can_sps_list[match(sps_order, can_sps_list)]
can_sps_list = can_sps_list[!is.na(can_sps_list)]


# retain good alignments for genes with annotated TE domains
ali = ali [ ali$evalue < 1e-5, ]
ali = merge(can, ali, by.x = "gene", by.y = "qseqid", all.x = T, all.y = F)
ali = ali[!duplicated(ali$gene),]

# merge
alc = merge(x=ali, y=tec, by.x="sseqid", by.y="TE_id", all.x=T, all.y = F)
alc = unique(alc[,c("gene","species","TE_class","architecture")])


#### General loop per species ####
gen_fusions = data.frame()

for (sps in can_sps_list) {
  
  print(sps)
  
  #### Define per-species candidates ####
  
  # create dataframe to store per-gene diagnostics
  spi_fusions = data.frame(gene = can[can$species == sps,"gene"], architecture = can[can$species == sps,"architecture"])
  spi_fusions$core_class_string = lapply(spi_fusions$gene, function(x) paste(sort(gec[gec$gene == x,"class"]), collapse = ","))
  spi_fusions$TE_class_string = lapply(spi_fusions$gene, function(x) paste(alc[alc$gene == x,"TE_class"], collapse = "|"))
  
  if (sps %in% val_sps_list) {
    
    # read pfam domain coordinates
    spi_pfam = read.table(sprintf("results_TEfusions/data/fusions.%s.Pfamscan.seqs.csv", sps), stringsAsFactors = F,
                          header = F, col.names = c("gene","dom_start","dom_end","pfam_id","pfam_dom","pfam_seq"))
    spi_pfam$is_transposon = spi_pfam$pfam_dom %in% ted$domain
    spi_pfam$is_core = spi_pfam$pfam_dom %in% cla_core_domains
    
    
    #### Assembly contiguity ####
    
    # read assembly contiguity
    # transcripts that are contiguous in the genome assembly (no NNN stretches within gene model)
    if (file.info(sprintf("results_TEfusions/data/fusions.%s.nonNstretch.txt", sps))$size > 0) {
      spi_acon = read.table(sprintf("results_TEfusions/data/fusions.%s.nonNstretch.txt", sps), stringsAsFactors = F,
                            header = F)[,1]
      spi_fusions$is_contiguous_assembly = spi_fusions$gene %in% spi_acon
    } else {
      spi_fusions$is_contiguous_assembly = F
    }
    
    
    #### Expression ####
    
    # read expression
    spi_expr = read.table(sprintf(
      "results_TEfusions/data/fusions.%s.salmon_quant.csv", sps), 
      stringsAsFactors = F, header = T)
    # read expression contiguity
    spi_exco = read.table(sprintf(
      "results_TEfusions/data/fusions.%s.salmon_cover.csv", sps),
      header = T, stringsAsFactors = F)
    spi_exco = merge(spi_exco, spi_expr, by.x="Transcript", by.y="Transcript")
    
    if (nrow(spi_expr) > 1 & nrow(spi_exco) >1) { 
      spi_expr_s = spi_expr[c("Gene","TPM","Sample")]
      spi_expr_m = reshape(spi_expr_s, idvar = "Gene", timevar = "Sample", direction = "wide")
      rownames(spi_expr_m) = spi_expr_m$Gene
      spi_expr_m = within(spi_expr_m, rm("Gene"))
      colnames(spi_expr_m) = stringr::str_split(colnames(spi_expr_m), pattern = "\\.", simplify = T)[,2]
      # number of samples in which each gene is expressed
      spi_expr_count = data.frame(gene = names(rowSums(spi_expr_m>0)), num_expressed=rowSums(spi_expr_m>0))
      spi_fusions = merge(spi_fusions, spi_expr_count, all.x = T, all.y = F)
      spi_fusions$num_expressed [is.na(spi_fusions$num_expressed )] = 0
      spi_fusions$is_expressed = spi_fusions$num_expressed>0
      spi_exco_gr = GenomicRanges::makeGRangesFromDataFrame(spi_exco, start.field = "start_pep", end.field = "end_pep", seqnames.field = "Gene")
      spi_full_coverage = c() # empty vector (gets filled below)
      
      
      #### Contiguous expression coverage ####
      
      # pfam coordinates per gene to identify core and TE region
      #### FIRST:
      # - core region: first core domain found in the gene
      # - TE region: last TE domain found in the gene
      # - result: the gap will be one of the most extensive in that gene that bridges the fusion
      spi_pfam_gene = data.frame(row.names = spi_fusions$gene)
      spi_pfam_gene$min_core = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        min(spi_pfam[spi_pfam$gene == x & spi_pfam$is_core, "dom_start" ])))
      spi_pfam_gene$max_core = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        min(spi_pfam[spi_pfam$gene == x & spi_pfam$is_core, "dom_end" ]) ))
      spi_pfam_gene$min_te = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        max(spi_pfam[spi_pfam$gene == x & spi_pfam$is_transposon, "dom_start" ]) ))
      spi_pfam_gene$max_te = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        max(spi_pfam[spi_pfam$gene == x & spi_pfam$is_transposon, "dom_end" ]) ))
      spi_pfam_gene = spi_pfam_gene[is.finite(rowSums(spi_pfam_gene)),]
      spi_pfam_gene$gene = rownames(spi_pfam_gene)
      
      # define GRanges of core region
      spi_pcor_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_pfam_gene,start.field = "min_core",end.field = "max_core", seqnames.field = "gene")
      # define GRanges of TE region
      spi_ptes_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_pfam_gene,start.field = "min_te",end.field = "max_te", seqnames.field = "gene")
      # define GRanges gaps between core and TE
      spi_pgap_gr = IRanges::pgap(spi_pcor_gr, spi_ptes_gr)
      
      # overlap between domains and coverage deserts
      spi_ovp_pgap_exco = IRanges::findOverlaps(query = spi_pgap_gr, spi_exco_gr, type = "within")
      spi_full_coverage = c(spi_full_coverage,names(spi_pgap_gr)[unique(spi_ovp_pgap_exco@from)])
      
      #### SECOND:
      # - core region: last core domain found in the gene
      # - TE region: first TE domain found in the gene
      # - result: the gap will be one of the most extensive in that gene that bridges the fusion
      spi_pfam_gene = data.frame(row.names = spi_fusions$gene)
      spi_pfam_gene$min_core = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        max(spi_pfam[spi_pfam$gene == x & spi_pfam$is_core, "dom_start" ])))
      spi_pfam_gene$max_core = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        max(spi_pfam[spi_pfam$gene == x & spi_pfam$is_core, "dom_end" ]) ))
      spi_pfam_gene$min_te = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        min(spi_pfam[spi_pfam$gene == x & spi_pfam$is_transposon, "dom_start" ]) ))
      spi_pfam_gene$max_te = unlist(lapply(rownames(spi_pfam_gene), function(x) 
        min(spi_pfam[spi_pfam$gene == x & spi_pfam$is_transposon, "dom_end" ]) ))
      spi_pfam_gene = spi_pfam_gene[is.finite(rowSums(spi_pfam_gene)),]
      spi_pfam_gene$gene = rownames(spi_pfam_gene)
      
      # define GRanges of core region
      spi_pcor_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_pfam_gene,start.field = "min_core",end.field = "max_core", seqnames.field = "gene")
      # define GRanges of TE region
      spi_ptes_gr = GenomicRanges::makeGRangesFromDataFrame(df = spi_pfam_gene,start.field = "min_te",end.field = "max_te", seqnames.field = "gene")
      # define GRanges gaps between core and TE
      spi_pgap_gr = IRanges::pgap(spi_pcor_gr, spi_ptes_gr)
      
      # overlap between domains and coverage deserts
      spi_ovp_pgap_exco = IRanges::findOverlaps(query = spi_pgap_gr, spi_exco_gr, type = "within")
      spi_full_coverage = c(spi_full_coverage,names(spi_pgap_gr)[unique(spi_ovp_pgap_exco@from)])
      
      # report full coverage
      spi_fusions$is_contiguous_coverage = spi_fusions$gene %in% spi_full_coverage
      
    } else {
      
      # if there are no expressed genes, record expression data as FALSE
      spi_fusions$num_expressed = NA
      spi_fusions$is_expressed = F
      spi_fusions$is_contiguous_coverage = F
      
    }
    
    # map gene to transcript names
    spi_ge2t = spi_expr
    spi_ge2t = unique(within(spi_ge2t, rm("Sample","Length","EffectiveLength","TPM","NumReads")))
    
    #### Gene structure ####
    
    spi_gest = read.table(sprintf(
      "results_TEfusions/data/fusions.%s.exons_per_gene.csv", sps), 
      stringsAsFactors = F, header = F, col.names = c("Transcript","num_exons"))
    spi_gest = merge(spi_gest, spi_ge2t, by.x="Transcript", by.y="Transcript", all.y=F, all.x=T)
    spi_gest_monoexonic = spi_gest[spi_gest$num_exons == 1, "Gene" ]
    spi_fusions$is_monoexonic = spi_fusions$gene %in% spi_gest_monoexonic
    
    #### Genome clusters ####
    # defined using bedtools cluster tool, at a fixed distance (-d 100000)
    
    spi_clus = read.table(sprintf(
      "results_TEfusions/data/fusions.%s.diamond.annot.bed",sps),
      stringsAsFactors = F, header = F, col.names = c("chr","start","end","Transcript","cluster"))
    spi_clus = merge(spi_clus, spi_ge2t, by.x="Transcript", by.y="Transcript", all.y=F, all.x=T)
    spi_cluster_counts = table(spi_clus$cluster)
    spi_cluster_list = names(spi_cluster_counts)[spi_cluster_counts>1]
    spi_cluster_genes = spi_clus[spi_clus$cluster %in% spi_cluster_list,"Gene"]
    spi_fusions$is_in_cluster = spi_fusions$gene %in% spi_cluster_genes
    
    # if species is in the validation list, record
    spi_fusions$is_validated = T
    
  } else {
    
    # if species is not in the validation list, record all as FALSE
    spi_fusions$is_contiguous_assembly = F
    spi_fusions$num_expressed = NA
    spi_fusions$is_expressed = F
    spi_fusions$is_contiguous_coverage = F
    spi_fusions$is_monoexonic = F
    spi_fusions$is_in_cluster = F
    
    # if species is not in the validation list, record
    spi_fusions$is_validated = F
    
  } 
  
  #### Store data ####
  gen_fusions = rbind(gen_fusions, spi_fusions)
  
}


#### Summary ####

# species order
gen_fusions$species = stringr::str_split(gen_fusions$gene, "_", simplify = T)[,1]
gen_fusions$species = factor(gen_fusions$species, levels = can_sps_list)

# high-level TE taxonomy
gen_fusions$TE_class_string = unlist(gen_fusions$TE_class_string)
gen_fusions$TE_class = "Unknown"
gen_fusions[grepl(x = gen_fusions$TE_class_string, pattern = "Class_II_DNA_Transposition"),"TE_class"] = "DNA"
gen_fusions[grepl(x = gen_fusions$TE_class_string, pattern = "Class_I_Retrotransposition"),"TE_class"] = "RT"
gen_fusions[grepl(x = gen_fusions$TE_class_string, pattern = "Class_I_Retrotransposition") & grepl(x = gen_fusions$TE_class_string, pattern = "Class_II_DNA_Transposition"),"TE_class"] = "DNA/RT"
# low-level TE taxonomy
gen_fusions$TE_type = unlist(lapply(stringr::str_split(gen_fusions$TE_class_string, ";"), function(x) {
  len = max(1, length(x)-1 )
  x[len]
}))
gen_fusions$TE_classtype = paste(gen_fusions$TE_class, gen_fusions$TE_type)
# if TE class or type are unknown, mix them
gen_fusions$TE_classtype[gen_fusions$TE_classtype == "Unknown NA" | gen_fusions$TE_classtype == "Unknown Interspersed_Repeat"] = "Unknown"

# TE domains
gen_fusions$TE_domain_string = unlist(
  lapply(stringr::str_split(gen_fusions$architecture, " "),
         function(x) paste(sort(unique(x[x %in% ted$domain])), collapse = ",")))
# core domains
gen_fusions$core_class_string = unlist(gen_fusions$core_class_string)
# simplified architecture
gen_fusions$architecture_simple = unlist(
  lapply(stringr::str_split(gen_fusions$architecture, " "),
         function(x) paste(unique(x),collapse = ",")))

# summarise expression status
gen_fusions$evidence_status = "no evidence"
gen_fusions[gen_fusions$is_expressed  & gen_fusions$is_contiguous_coverage  & gen_fusions$is_contiguous_assembly  & gen_fusions$is_validated,"evidence_status"] = "contiguous assembly, complete expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_coverage & gen_fusions$is_contiguous_assembly  & gen_fusions$is_validated,"evidence_status"] = "contiguous assembly, partial expression"
gen_fusions[!gen_fusions$is_expressed & !gen_fusions$is_contiguous_coverage & gen_fusions$is_contiguous_assembly  & gen_fusions$is_validated,"evidence_status"] = "contiguous assembly, no expression"
gen_fusions[!gen_fusions$is_expressed & gen_fusions$is_contiguous_coverage  & gen_fusions$is_contiguous_assembly  & gen_fusions$is_validated,"evidence_status"] = "contiguous assembly, no expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_assembly & gen_fusions$is_contiguous_assembly  & gen_fusions$is_validated,"evidence_status"] = "broken assembly, complete expression"
gen_fusions[gen_fusions$is_expressed  & !gen_fusions$is_contiguous_assembly & !gen_fusions$is_contiguous_assembly & gen_fusions$is_validated,"evidence_status"] = "broken assembly, partial expression"
gen_fusions[!gen_fusions$is_expressed & !gen_fusions$is_contiguous_assembly & gen_fusions$is_validated,"evidence_status"] = "broken assembly, no expression"
gen_fusions[gen_fusions$species %in% sps_transcriptome_only, "evidence_status"] = "contiguous assembly, complete expression"


# factors
gen_fusions$evidence_status = factor(
  gen_fusions$evidence_status,
  levels = c("contiguous assembly, complete expression",
             "contiguous assembly, partial expression",
             "contiguous assembly, no expression",
             "broken assembly, complete expression",
             "broken assembly, partial expression",
             "broken assembly, no expression",
             "no evidence"))

# summarise genome structure status
gen_fusions$structure_status = "unknown"
gen_fusions[gen_fusions$is_monoexonic & gen_fusions$is_in_cluster,"structure_status"]   = "monoexon, cluster"
gen_fusions[gen_fusions$is_monoexonic & !gen_fusions$is_in_cluster,"structure_status"]  = "monoexon, no cluster"
gen_fusions[!gen_fusions$is_monoexonic & gen_fusions$is_in_cluster,"structure_status"]  = "multiexon, cluster"
gen_fusions[!gen_fusions$is_monoexonic & !gen_fusions$is_in_cluster,"structure_status"] = "multiexon, no cluster"
gen_fusions$structure_status = factor(
  gen_fusions$structure_status,
  levels=c("monoexon, cluster","monoexon, no cluster","multiexon, cluster","multiexon, no cluster","unknown"))

# simplify core class string (mix Chromo+SNF+PHD into single category)
gen_fusions$core_class_string_simple = gen_fusions$core_class_string
gen_fusions$core_class_string_simple[grepl("Chromo,", gen_fusions$core_class_string_simple)] = "Chromo+SNF2_N+(PHD+others)"
simple_string_others = names(which(table(gen_fusions$core_class_string_simple) < 5))
gen_fusions$core_class_string_simple[gen_fusions$core_class_string_simple %in% simple_string_others] = "Others"


# barplots with basic counts
pdf("results_TEfusions/gene_fusion_evidence-barplot.pdf", height = 8, width = 6)
par(mar = c(5.1, 20, 4.1, 2.1))
barplot(
  table(gen_fusions$evidence_status), horiz = T, las=1,
  names.arg = paste(names(table(gen_fusions$evidence_status)),"n =", table(gen_fusions$evidence_status)),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)

barplot(
  table(gen_fusions$species), horiz = T, las=1,
  names.arg = paste(names(table(gen_fusions$species)),"n =", table(gen_fusions$species)),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.3)

barplot(
  table(gen_fusions$TE_class), horiz = T, las=1,
  names.arg = paste(names(table(gen_fusions$TE_class)),"n =", table(gen_fusions$TE_class)),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)

barplot(
  table(gen_fusions$TE_classtype), horiz = T, las=1,
  names.arg = paste(names(table(gen_fusions$TE_classtype)),"n =", table(gen_fusions$TE_classtype)),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)

barplot(
  table(gen_fusions$core_class_string), horiz = T, las=1,
  names.arg = paste(names(table(gen_fusions$core_class_string)),"n =", table(gen_fusions$core_class_string)),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)


# order of simple string cross-tabulation
class_string_sorted_tab = sort(table(gen_fusions$core_class_string_simple))
class_string_sorted_nam = names(class_string_sorted_tab)
barplot(
  class_string_sorted_tab, horiz = T, las=1,
  names.arg = paste(names(class_string_sorted_tab),"n =", class_string_sorted_tab),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)

# how many hits to each type of transposon?
xtab = t(table(gen_fusions$core_class_string_simple, gen_fusions$TE_class)[class_string_sorted_nam,])
barplot(
  xtab, 
  horiz = T, las=1, col=c("red","orange","gray"),
  names.arg = paste(names(class_string_sorted_tab)," n=", xtab[1,], "/",xtab[2,], "/",xtab[3,], sep=""),
  main = sprintf("TE class similarity\nn=%i", nrow(gen_fusions)),cex.names = 0.8)
legend("bottomright", fill=c("red","orange","gray"), legend = c("DNA","RT","Unk"), bty = "n")

# how many species per combination of core domains?
xtab = table(gen_fusions$core_class_string_simple, gen_fusions$species)[class_string_sorted_nam,]
xtab = apply(xtab, 1, function(x) sum(x>0))
barplot(
  xtab, 
  horiz = T, las=1, 
  names.arg = paste(names(xtab),"n =", xtab),
  main = sprintf("num species having each fusion core combination\nn=%i", sum(xtab)),cex.names = 0.8)

# # barplots segmented by core domain
list_gene_fam = cla$gene_fam
gen_fusions_counts_per_core = unlist(lapply(list_gene_fam, function(x) length(grep(x, gen_fusions$core_class_string))))
names(gen_fusions_counts_per_core) = list_gene_fam
gen_fusions_counts_per_core = gen_fusions_counts_per_core[gen_fusions_counts_per_core>0]
barplot(
  gen_fusions_counts_per_core, horiz = T, las=1,
  names.arg = paste(names(gen_fusions_counts_per_core),"n =", gen_fusions_counts_per_core),
  main = sprintf("TE fusions\nn=%i", nrow(gen_fusions)),cex.names = 0.8)

dev.off()


# cross-tablulations of evidence status and fusion class
pdf("results_TEfusions/gene_fusion_evidence-table.pdf", height = 3, width = 3)

col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))
xtab = table(gen_fusions$core_class_string_simple, gen_fusions$evidence_status)[rev(class_string_sorted_nam),]
pheatmap(xtab,  
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         main = sprintf("evidence sum per class|n=%i", sum(xtab)),
         border_color = "white", display_numbers = T, number_format = "%i")

xtab = table(gen_fusions$core_class_string_simple, gen_fusions$structure_status)[rev(class_string_sorted_nam),]
pheatmap(xtab,  
         color = col_blue(10), breaks = seq(0,10,length.out = 11)-0.01, 
         cellwidth = 4, cellheight = 4, na_col = "grey", 
         cluster_rows = F, cluster_cols = F, 
         fontsize = 5, 
         main = sprintf("structure sum per class|n=%i", sum(xtab)),
         border_color = "white", display_numbers = T, number_format = "%i")


dev.off()


# write output
write.table(
  gen_fusions, 
  file = "results_TEfusions/gene_fusion_evidence.csv" ,
  quote = F, sep = "\t", row.names = F)


#### Add TE orthogroup ####

# plot pairwise orthologies of gene and TE domain
pdf("results_TEfusions/retrotransposon_OGs_overlap.pdf", height = 10, width = 6)
for (cli in list_gene_fam) {
  
  
  print(cli)
  # subset
  gen_fusions_i = gen_fusions[grepl(gen_fusions$core_class_string, pattern = cli) ,]
  # remove rows that contain Chromo domains in case we're not analysing the Chromo domain
  if (cli != "Chromo") {
    gen_fusions_i = gen_fusions_i[!grepl(gen_fusions_i$core_class_string, pattern = "Chromo"),]
  }
  gen_fusions_i_order = gen_fusions_i$gene
  
  
  # add orthology group (only within said gene family)
  ort_i = ort[ort$family == cli,c("gene","orthogroup")]
  gen_fusions_i = merge(gen_fusions_i, ort_i, by.x = "gene", by.y="gene", all.x=T, all.y=F)
  # gen_fusions_i$orthogroup[is.na(gen_fusions_i$orthogroup)] = "NA"
  gen_fusions_i$orthogroup = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
  gen_fusions_i = gen_fusions_i[match(gen_fusions_i_order, gen_fusions_i$gene),]
  
  # order rows according to combination of orthogroup (first) and species (second)
  gen_fusions_i$species = droplevels(gen_fusions_i$species)
  gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$orthogroup)
  gen_fusions_i = gen_fusions_i[order(gen_fusions_i$orthogroup, gen_fusions_i$species),]
  species_core_order = unique(gen_fusions_i$species_core)
  
  # add TE orthology
  gen_fusions_i = merge(gen_fusions_i, orte_sum, by.x="gene", by.y="gene", all.x=T, all.y=F)
  # gen_fusions_i[is.na(gen_fusions_i$orthogroup_TE),"orthogroup_TE"] = "NA"
  
  # remove DNA or unknown TEs
  # gen_fusions_i = gen_fusions_i[gen_fusions_i$TE_class == "RT",]
  gen_fusions_i = gen_fusions_i[!is.na(gen_fusions_i$orthogroup_TE) & !is.na(gen_fusions_i$orthogroup),]
  
  if (nrow(gen_fusions_i)>0) {
    
    # table TE classification
    xtab = as.data.frame(table(OG = gen_fusions_i$orthogroup,TEG = gen_fusions_i$orthogroup_TE))
    xtab = xtab[xtab$Freq>0,]
    if(nrow(xtab)>1){
      alluvial(
        xtab[,1:2], 
        freq = xtab[,3],
        col = c("lightblue"), gap.width = 0.1, blocks = T, 
        border=NA, cex=0.6)
    } else {
      plot(0,0, xlab = "", ylab = "", frame.plot = F, axes = F, col=NA, pch=19)
      text(0,0,sprintf("%s is %s", as.character(xtab$OG[1]), as.character(xtab$TEG[1])))
    }
    title(main=sprintf("%s RT domain fusions",cli),
          sub=sprintf("n=%i", nrow(gen_fusions_i)))
    
  }
  
}
dev.off()

pdf("results_TEfusions/retrotransposon_OGs_overlap_detail.pdf", height = 10, width = 6)
for (cli in list_gene_fam) {
  
  
  print(cli)
  # subset
  gen_fusions_i = gen_fusions[grepl(gen_fusions$core_class_string, pattern = cli) ,]
  # remove rows that contain Chromo domains in case we're not analysing the Chromo domain
  if (cli != "Chromo") {
    gen_fusions_i = gen_fusions_i[!grepl(gen_fusions_i$core_class_string, pattern = "Chromo"),]
  }
  gen_fusions_i_order = gen_fusions_i$gene
  
  
  # add orthology group (only within said gene family)
  ort_i = ort[ort$family == cli,c("gene","orthogroup")]
  gen_fusions_i = merge(gen_fusions_i, ort_i, by.x = "gene", by.y="gene", all.x=T, all.y=F)
  # gen_fusions_i$orthogroup[is.na(gen_fusions_i$orthogroup)] = "NA"
  gen_fusions_i$orthogroup = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
  gen_fusions_i = gen_fusions_i[match(gen_fusions_i_order, gen_fusions_i$gene),]
  
  # order rows according to combination of orthogroup (first) and species (second)
  gen_fusions_i$species = droplevels(gen_fusions_i$species)
  gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$orthogroup)
  gen_fusions_i = gen_fusions_i[order(gen_fusions_i$orthogroup, gen_fusions_i$species),]
  species_core_order = unique(gen_fusions_i$species_core)
  
  # add TE orthology
  gen_fusions_i = merge(gen_fusions_i, orte_sum, by.x="gene", by.y="gene", all.x=T, all.y=F)
  # gen_fusions_i[is.na(gen_fusions_i$orthogroup_TE),"orthogroup_TE"] = "NA"
  
  # remove DNA or unknown TEs
  # gen_fusions_i = gen_fusions_i[gen_fusions_i$TE_class == "RT",]
  gen_fusions_i = gen_fusions_i[!is.na(gen_fusions_i$orthogroup_TE) & !is.na(gen_fusions_i$orthogroup),]
  
  ogi_list = unique(gen_fusions_i$orthogroup)
  
  for (ogi in ogi_list) {
    
    gen_fusions_i_ogi = gen_fusions_i[gen_fusions_i$orthogroup == ogi,]
    if (nrow(gen_fusions_i_ogi)>0) {
      
      # table TE classification
      xtab = as.data.frame(table(OG = gen_fusions_i_ogi$species_core,TEG = gen_fusions_i_ogi$orthogroup_TE))
      xtab = xtab[xtab$Freq>0,]
      if(nrow(xtab)>1){
        alluvial(
          xtab[,1:2], 
          freq = xtab[,3],
          col = c("lightblue"), gap.width = 0.1, blocks = T, 
          border=NA, cex=0.6)
      } else {
        plot(0,0, xlab = "", ylab = "", frame.plot = F, axes = F, col=NA, pch=19)
        text(0,0,sprintf("%s is %s", as.character(xtab$OG[1]), as.character(xtab$TEG[1])))
      }
      title(main=sprintf("%s RT domain fusions\n%s",cli,ogi),
            sub=sprintf("n=%i", nrow(gen_fusions_i_ogi)))
      
    }
  }
  
}
dev.off()





#### Plot Trees with TE orthology ####

tree_fo = "gene_trees/"
for (cli in list_gene_fam) {
  
  print(cli)
  # subset
  gen_fusions_i = gen_fusions[grepl(gen_fusions$core_class_string, pattern = cli) ,]
  # remove rows that contain Chromo domains in case we're not analysing the Chromo domain
  if (cli != "Chromo") {
    gen_fusions_i = gen_fusions_i[!grepl(gen_fusions_i$core_class_string, pattern = "Chromo"),]
  }
  gen_fusions_i_order = gen_fusions_i$gene
  
  
  # add orthology group (only within said gene family)
  ort_i = ort[ort$family == cli,c("gene","orthogroup")]
  gen_fusions_i = merge(gen_fusions_i, ort_i, by.x = "gene", by.y="gene", all.x=T, all.y=F)
  
  # add TE orthology
  gen_fusions_i = merge(gen_fusions_i, orte_sum, by.x="gene", by.y="gene", all.x=T, all.y=F)
  
  # remove DNA or unknown TEs
  # gen_fusions_i = gen_fusions_i[gen_fusions_i$TE_class == "RT",]
  gen_fusions_i = gen_fusions_i[!is.na(gen_fusions_i$orthogroup_TE) & !is.na(gen_fusions_i$orthogroup),]
  
  if (nrow(gen_fusions_i)>0){
    
    # retrieve list of homology groups
    gen_fusions_i$homology_group_d = stringr::str_split(gen_fusions_i$orthogroup, pattern = "\\.", simplify = T)[,1]
    gen_fusions_i$homology_group_h = stringr::str_split(gen_fusions_i$orthogroup, pattern = "\\.", simplify = T)[,2]
    gen_fusions_i$homology_group = paste(gen_fusions_i$homology_group_d, gen_fusions_i$homology_group_h, sep=".")
    # gen_fusions_i$orthogroup_id = stringr::str_split(gen_fusions_i$orthogroup, pattern = ":", simplify = T)[,1]
    gen_fusions_i$orthogroup_id = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
    
    # FIRST
    # plot core trees with TE domain orthology (RVT and rve domains)
    hgs_w_TE = unique(gen_fusions_i$homology_group)
    pdf(file = sprintf("results_TEfusions/tree_core_with_RT_domains.%s.pdf", cli),width = 12, height = 12)
    for (ogi in hgs_w_TE) {
      
      if (file.exists(sprintf("%s/euk.%s.seqs.iqtree.treefile", tree_fo, ogi))) {
        
        gphy = read.tree(
          file=sprintf("%s/euk.%s.seqs.iqtree.treefile", tree_fo, ogi))
        
        # shorten tip labels
        gphy$tip.label = stringr::str_split(gphy$tip.label, pattern = "\\|", simplify = T)[,1]
        
        # classes per fusion
        gphy_fus = merge(gen_fusions_i, data.frame(tip=gphy$tip.label), all.x=F, all.y=T, by.x="gene", by.y="tip")
        rownames(gphy_fus) = gphy_fus$gene
        gphy_fus = gphy_fus[gphy$tip.label,]
        gphy_fus_tab = table(gphy_fus$orthogroup_TE)
        gphy_fus_cla = names(gphy_fus_tab)
        gphy_fus_col = rainbow(n=length(gphy_fus_cla), end = 0.8, v = 0.8)
        gphy_fus$orthogroup_TE[is.na(gphy_fus$orthogroup_TE) & !is.na(gphy_fus$core_class_string)] = "other"
        gphy_fus$orthogroup_TE[is.na(gphy_fus$orthogroup_TE) & is.na(gphy_fus$core_class_string)] = "none"
        
        # plot gene tree      
        text_size=0.3
        gphy$edge.length[gphy$edge.length<0.01] = 0.01
        plot.phylo(gphy, font=1, type="u", label.offset = 0.01, cex=text_size, edge.color = "gray",
                   root.edge = T, show.tip.label = F, align.tip.label = F, underscore = T)
        # label tips        
        unlabeled_tips = which(gphy_fus$orthogroup_TE == "none")
        tiplabels(text=gphy$tip.label[unlabeled_tips], tip = unlabeled_tips, col = alpha("darkgray",0.7), bg = NA, frame = "none", cex=0.4)
        # plot classes per fusion
        for (fui in seq_along(gphy_fus_cla)) {
          gphy_fus_cla_i = gphy_fus_cla[fui]
          gphy_fus_col_i = gphy_fus_col[fui]
          tips_w_TE = which(gphy_fus$orthogroup_TE == gphy_fus_cla_i)
          tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.7, tip = tips_w_TE)
          tiplabels(text=gphy$tip.label[tips_w_TE], tip = tips_w_TE, col = alpha(gphy_fus_col_i,0.8), bg = NA, frame = "none", cex=0.5)
        }
        # mark tips with full evidence
        tips_w_evidence = which(gphy_fus$evidence_status == "contiguous assembly, complete expression")
        tiplabels(pch = 1, col = "black", height = 4, cex = 1, tip = tips_w_evidence)
        legend("bottomright", col = gphy_fus_col, legend = paste(names(gphy_fus_tab),"n =", gphy_fus_tab), bty="n", pch=19, cex=0.7)
        title(
          main = sprintf("%s",ogi),
          sub =  sprintf("%i genes\nn=%i fusion(s) with RT domains\n(%i pass filters)", length(gphy$tip.label), sum(gphy_fus$orthogroup_TE!="none"), length(tips_w_evidence)))
        # add scale bar
        add.scale.bar(x=0, y=-5)
        
      }
      
    }
    dev.off()
    
    
    
    # SECOND
    # plot TE trees with core orthology (RVT and rve domains)
    TEs_w_OG_d = str_split((unlist(str_split(gen_fusions_i$orthogroup_TE,","))), "\\.", simplify = T)[,1]
    TEs_w_OG_h = str_split((unlist(str_split(gen_fusions_i$orthogroup_TE,","))), "\\.", simplify = T)[,2]
    TEs_w_OG = sort(unique(paste(TEs_w_OG_d, TEs_w_OG_h, sep = ".")))
    pdf(file = sprintf("results_TEfusions/tree_RT_domains_with_core.%s.pdf", cli),width = 12, height = 12)
    for (ogi in TEs_w_OG) {
      
      if (file.exists(sprintf("%s/euk.%s.seqs.iqtree.treefile", orte_fn, ogi))) {
        
        gphy = read.tree(
          file=sprintf("%s/euk.%s.seqs.iqtree.treefile", orte_fn, ogi))
        
        # shorten tip labels
        gphy$tip.label = stringr::str_split(gphy$tip.label, pattern = "\\|", simplify = T)[,1]
        # remove domain id
        gphy$tip.label = gsub(pattern = "_\\d+-\\d+$", replacement="", gphy$tip.label)
        
        # classes per fusion
        gphy_fus = merge(gen_fusions_i, data.frame(tip=gphy$tip.label), all.x=F, all.y=T, by.x="gene", by.y="tip")
        gphy_fus = gphy_fus[match(gphy$tip.label, gphy_fus$gene),]
        gphy_fus_tab = table(gphy_fus$orthogroup_id)
        gphy_fus_cla = names(gphy_fus_tab)
        gphy_fus_col = rainbow(n=length(gphy_fus_cla), end = 0.8, v = 0.8)
        gphy_fus$orthogroup_id[is.na(gphy_fus$orthogroup_id)] = "none"
        
        # plot gene tree      
        text_size=0.3
        gphy$edge.length[gphy$edge.length<0.01] = 0.01
        plot.phylo(gphy, font=1, type="u", label.offset = 0.01, cex=text_size, edge.color = "gray",
                   root.edge = T, show.tip.label = T, align.tip.label = F, underscore = T,tip.color = alpha("gray10",0.8))
        # plot classes per fusion
        for (fui in seq_along(gphy_fus_cla)) {
          gphy_fus_cla_i = gphy_fus_cla[fui]
          gphy_fus_col_i = gphy_fus_col[fui]
          tips_w_TE = which(gphy_fus$orthogroup_id == gphy_fus_cla_i)
          tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.8, tip = tips_w_TE)
        }
        tips_w_evidence = which(gphy_fus$evidence_status == "contiguous assembly, complete expression")
        tiplabels(pch = 1, col = "black", height = 4, cex = 1, tip = tips_w_evidence)
        legend("bottomright", col = gphy_fus_col, legend = paste(names(gphy_fus_tab),"n =", gphy_fus_tab), bty="n", pch=19, cex=0.7)
        title(
          main = sprintf("%s",ogi),
          sub =  sprintf("%i genes\nn=%i fusion(s) with %s domains\n(%i pass filters)", length(gphy$tip.label), sum(gphy_fus$orthogroup_id != "none"), cli, length(tips_w_evidence)))
        # add scale bar
        add.scale.bar(x=0, y=-5)
        
      }
      
    }
    dev.off()
    
  }
  
}







#### Heatmaps counts ####

pdf("results_TEfusions/heatmap_TE_fusion_evidence.pdf", height = 9, width = 6)
for (cli in list_gene_fam) {
  
  # subset    
  gen_fusions_i = gen_fusions[grepl(gen_fusions$core_class_string, pattern = cli) ,]
  # remove rows that contain Chromo domains in case we're not analysing the Chromo domain
  if (cli != "Chromo") {
    gen_fusions_i = gen_fusions_i[!grepl(gen_fusions_i$core_class_string, pattern = "Chromo"),]
  }
  gen_fusions_i_order = gen_fusions_i$gene
  
  # add orthology group (only within said gene family)
  ort_i = ort[ort$family == cli,c("gene","orthogroup")]
  gen_fusions_i = merge(gen_fusions_i, ort_i, by.x = "gene", by.y="gene", all.x=T, all.y=F)
  gen_fusions_i$orthogroup[is.na(gen_fusions_i$orthogroup)] = "NA"
  # gen_fusions_i$orthogroup = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
  gen_fusions_i = gen_fusions_i[match(gen_fusions_i_order, gen_fusions_i$gene),]
  
  # # order rows according to combinations of species and core domains (proxy to gene family)
  # gen_fusions_i$species = droplevels(gen_fusions_i$species)
  # gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$core_class_string)
  # species_core_order = unique(gen_fusions_i$species_core)
  
  # order rows according to combination of orthogroup (first) and species (second)
  gen_fusions_i$species = droplevels(gen_fusions_i$species)
  gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$orthogroup)
  gen_fusions_i = gen_fusions_i[order(gen_fusions_i$orthogroup, gen_fusions_i$species),]
  species_core_order = unique(gen_fusions_i$species_core)
  
  # cross-tabulations...    
  if (nrow(gen_fusions_i)>0) {
    # col separation
    colsep = c()
    
    # table TE classification
    tab = xtabs(formula = ~ species_core + TE_classtype ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = tam
    colsep = c(colsep, ncol(tat))
    
    # table TE domains
    gen_fusions_i_TEdoms_list = sort(unique(stringr::str_split(paste(unique(gen_fusions_i$TE_domain_string), collapse = ","), ",")[[1]]))
    for (dom in gen_fusions_i_TEdoms_list) {
      gen_fusions_j = gen_fusions_i[,c("species_core","TE_domain_string")]
      gen_fusions_j$domain = grepl(pattern = sprintf("\\b%s\\b", dom), gen_fusions_j$TE_domain_string) * 1
      tab = aggregate(gen_fusions_j$domain, by=list(gen_fusions_j$species_core), sum)
      # rownames(tab) = tab$Group.1
      # tab = tab[species_core_order,]
      tat = cbind(tat, tab$x)
      colnames(tat)[length(colnames(tat))] = dom
    }
    colsep = c(colsep, ncol(tat))
    
    # table evidence status
    tab = xtabs(formula = ~ species_core + evidence_status ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = cbind(tat, tam)
    colsep = c(colsep, ncol(tat))
    
    # table monexon
    tab = xtabs(formula = ~ species_core + structure_status ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = cbind(tat, tam)
    colsep = c(colsep, ncol(tat))
    
    # plot      
    tas = matrix(tat[species_core_order,], nrow = length(species_core_order))
    rownames(tas) = species_core_order
    colnames(tas) = colnames(tat)
    pheatmap(tas,  
             gaps_col = colsep,
             color = col_blue(5), breaks = seq(0,5,length.out = 6)-0.01, 
             cellwidth = 4, cellheight = 4, na_col = "grey", 
             cluster_rows = F, cluster_cols = F, 
             labels_row = stringr::str_trunc(rownames(tas), width = 40),
             fontsize = 5, 
             main = sprintf("%s TEs", cli),
             border_color = "white", display_numbers = T, number_format = "%i")
  }
}
dev.off()

pdf("results_TEfusions/heatmap_TE_fusion_homology.pdf", height = 9, width = 6)
for (cli in list_gene_fam) {
  
  # subset    
  gen_fusions_i = gen_fusions[grepl(gen_fusions$core_class_string, pattern = cli) ,]
  # remove rows that contain Chromo domains in case we're not analysing the Chromo domain
  if (cli != "Chromo") {
    gen_fusions_i = gen_fusions_i[!grepl(gen_fusions_i$core_class_string, pattern = "Chromo"),]
  }
  gen_fusions_i_order = gen_fusions_i$gene
  
  # add orthology group (only within said gene family)
  ort_i = ort[ort$family == cli,c("gene","orthogroup")]
  gen_fusions_i = merge(gen_fusions_i, ort_i, by.x = "gene", by.y="gene", all.x=T, all.y=F)
  gen_fusions_i$orthogroup[is.na(gen_fusions_i$orthogroup)] = "NA"
  # gen_fusions_i$orthogroup = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
  gen_fusions_i = gen_fusions_i[match(gen_fusions_i_order, gen_fusions_i$gene),]
  
  # # order rows according to combinations of species and core domains (proxy to gene family)
  # gen_fusions_i$species = droplevels(gen_fusions_i$species)
  # gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$core_class_string)
  # species_core_order = unique(gen_fusions_i$species_core)
  
  # order rows according to combination of orthogroup (first) and species (second)
  gen_fusions_i$species = droplevels(gen_fusions_i$species)
  gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$orthogroup)
  gen_fusions_i = gen_fusions_i[order(gen_fusions_i$orthogroup, gen_fusions_i$species),]
  species_core_order = unique(gen_fusions_i$species_core)
  
  # cross-tabulations...    
  if (nrow(gen_fusions_i)>0) {
    # col separation
    colsep = c()
    
    # table TE classification
    tab = xtabs(formula = ~ species_core + TE_classtype ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = tam
    colsep = c(colsep, ncol(tat))
    
    # table TE domains
    gen_fusions_i_TEdoms_list = sort(unique(stringr::str_split(paste(unique(gen_fusions_i$TE_domain_string), collapse = ","), ",")[[1]]))
    for (dom in gen_fusions_i_TEdoms_list) {
      gen_fusions_j = gen_fusions_i[,c("species_core","TE_domain_string")]
      gen_fusions_j$domain = grepl(pattern = sprintf("\\b%s\\b", dom), gen_fusions_j$TE_domain_string) * 1
      tab = aggregate(gen_fusions_j$domain, by=list(gen_fusions_j$species_core), sum)
      # rownames(tab) = tab$Group.1
      # tab = tab[species_core_order,]
      tat = cbind(tat, tab$x)
      colnames(tat)[length(colnames(tat))] = dom
    }
    colsep = c(colsep, ncol(tat))
    
    # table evidence status
    tab = xtabs(formula = ~ species_core + evidence_status ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = cbind(tat, tam)
    colsep = c(colsep, ncol(tat))
    
    # table monexon
    tab = xtabs(formula = ~ species_core + structure_status ,
                data = gen_fusions_i, drop.unused.levels	= F)
    tam = matrix(tab, nrow = dim(tab)[1])
    rownames(tam) = rownames(tab)
    colnames(tam) = colnames(tab)
    tat = cbind(tat, tam)
    colsep = c(colsep, ncol(tat))
    
    # plot      
    tas = matrix(tat[species_core_order,], nrow = length(species_core_order))
    rownames(tas) = species_core_order
    colnames(tas) = colnames(tat)
    pheatmap(tas,  
             gaps_col = colsep,
             color = col_blue(5), breaks = seq(0,5,length.out = 6)-0.01, 
             cellwidth = 4, cellheight = 4, na_col = "grey", 
             cluster_rows = F, cluster_cols = F, 
             labels_row = stringr::str_trunc(rownames(tas), width = 40),
             fontsize = 5, 
             main = sprintf("%s TEs", cli),
             border_color = "white", display_numbers = T, number_format = "%i")
  }
}
dev.off()


# 
# # same but only combinations with full support
# 
# pdf("results_TEfusions/heatmap_TE_fusion_evidence_FULL.pdf", height = 10, width = 10)
# 
# # subset    
# gen_fusions_i = gen_fusions[gen_fusions$evidence_status == "contiguous assembly, complete expression" ,]
# gen_fusions_i_order = gen_fusions_i$gene
# 
# # add orthology group
# gen_fusions_i = merge(gen_fusions_i, ort, by.x = "gene", by.y="gene", all.x=T, all.y=F)
# gen_fusions_i$orthogroup[is.na(gen_fusions_i$orthogroup)] = "NA"
# # gen_fusions_i$orthogroup = stringr::str_trunc(gen_fusions_i$orthogroup, width = 40)
# gen_fusions_i = gen_fusions_i[match(gen_fusions_i_order, gen_fusions_i$gene),]
# 
# # order rows according to combination of orthogroup (first) and species (second)
# gen_fusions_i$species = droplevels(gen_fusions_i$species)
# gen_fusions_i$species_core = paste(gen_fusions_i$species, gen_fusions_i$orthogroup)
# gen_fusions_i = gen_fusions_i[order(gen_fusions_i$orthogroup, gen_fusions_i$species),]
# species_core_order = unique(gen_fusions_i$species_core)
# 
# # cross-tabulations...    
# if (nrow(gen_fusions_i)>0) {
#   # col separation
#   colsep = c()
#   
#   # table TE classification
#   tab = xtabs(formula = ~ species_core + TE_classtype ,
#               data = gen_fusions_i, drop.unused.levels	= F)
#   tam = matrix(tab, nrow = dim(tab)[1])
#   rownames(tam) = rownames(tab)
#   colnames(tam) = colnames(tab)
#   tat = tam
#   colsep = c(colsep, ncol(tat))
#   
#   # table TE domains
#   gen_fusions_i_TEdoms_list = sort(unique(stringr::str_split(paste(unique(gen_fusions_i$TE_domain_string), collapse = ","), ",")[[1]]))
#   for (dom in gen_fusions_i_TEdoms_list) {
#     gen_fusions_j = gen_fusions_i[,c("species_core","TE_domain_string")]
#     gen_fusions_j$domain = grepl(pattern = sprintf("\\b%s\\b", dom), gen_fusions_j$TE_domain_string) * 1
#     tab = aggregate(gen_fusions_j$domain, by=list(gen_fusions_j$species_core), sum)
#     # rownames(tab) = tab$Group.1
#     # tab = tab[species_core_order,]
#     tat = cbind(tat, tab$x)
#     colnames(tat)[length(colnames(tat))] = dom
#   }
#   colsep = c(colsep, ncol(tat))
#   
#   # table evidence status
#   tab = xtabs(formula = ~ species_core + evidence_status ,
#               data = gen_fusions_i, drop.unused.levels	= F)
#   tam = matrix(tab, nrow = dim(tab)[1])
#   rownames(tam) = rownames(tab)
#   colnames(tam) = colnames(tab)
#   tat = cbind(tat, tam)
#   colsep = c(colsep, ncol(tat))
#   
#   # table monexon
#   tab = xtabs(formula = ~ species_core + structure_status ,
#               data = gen_fusions_i, drop.unused.levels	= F)
#   tam = matrix(tab, nrow = dim(tab)[1])
#   rownames(tam) = rownames(tab)
#   colnames(tam) = colnames(tab)
#   tat = cbind(tat, tam)
#   colsep = c(colsep, ncol(tat))
#   
#   # plot      
#   tas = matrix(tat[species_core_order,], nrow = length(species_core_order))
#   rownames(tas) = species_core_order
#   colnames(tas) = colnames(tat)
#   pheatmap(tas,  
#            gaps_col = colsep,
#            color = col_blue(5), breaks = seq(0,5,length.out = 6)-0.01, 
#            cellwidth = 4, cellheight = 4, na_col = "grey", 
#            cluster_rows = F, cluster_cols = F, 
#            labels_row = stringr::str_trunc(rownames(tas), width = 100),
#            fontsize = 5, 
#            main = sprintf("All fusions (complete evidence only)", cli),
#            border_color = "white", display_numbers = T, number_format = "%i")
# }
# dev.off()
# 
# 
# print("Done!")
