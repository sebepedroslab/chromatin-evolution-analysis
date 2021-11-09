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


# read output
gen_fusions = read.table(
  file = "results_TEfusions/gene_fusion_evidence.csv" ,
  sep = "\t", header = T)



#### Plot Trees with TE orthology ####

tree_fo = "gene_trees/"
hgs_w_TE_list_d = str_split(list.files("*.treefile", path = "gene_trees_TEcore"), "\\.", simplify = T)[,2]
hgs_w_TE_list_h = str_split(list.files("*.treefile", path = "gene_trees_TEcore"), "\\.", simplify = T)[,3]
hgs_w_TE_list = paste(hgs_w_TE_list_d, hgs_w_TE_list_h, sep=".")

for (cli in cla$gene_fam) {
  
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
    hgs_w_TE = hgs_w_TE_list[hgs_w_TE_list_d == cli]
    if (length(hgs_w_TE)>0) {
      pdf(file = sprintf("results_TEfusions/tree_domain_with_RT_domains.%s.pdf", cli),width = 12, height = 12)
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
          gphy_fus$orthogroup_TE[is.na(gphy_fus$orthogroup_TE)] = "none"
          
          # plot gene tree      
          text_size=0.3
          gphy$edge.length[gphy$edge.length<0.01] = 0.01
          plot.phylo(gphy, font=1, type="u", label.offset = 0.01, cex=text_size, edge.color = "gray",
                     root.edge = T, show.tip.label = T, align.tip.label = F, underscore = T,tip.color = alpha("gray10",0.7),
                     main = sprintf("%s",ogi),sub=sprintf("n=%i fusion(s) with RT domains", sum(gphy_fus$orthogroup_TE!="none")))
          # plot classes per fusion
          for (fui in seq_along(gphy_fus_cla)) {
            gphy_fus_cla_i = gphy_fus_cla[fui]
            gphy_fus_col_i = gphy_fus_col[fui]
            tips_w_TE = which(gphy_fus$orthogroup_TE == gphy_fus_cla_i)
            tiplabels(pch = 19, col = gphy_fus_col_i, height = 4, cex = 0.7, tip = tips_w_TE)
          }
          tips_w_evidence = which(gphy_fus$evidence_status == "contiguous assembly, complete expression")
          tiplabels(pch = 1, col = "black", height = 4, cex = 1, tip = tips_w_evidence)
          legend("bottomright", col = gphy_fus_col, legend = paste(names(gphy_fus_tab),"n =", gphy_fus_tab), bty="n", pch=19, cex=0.7)
          # add scale bar
          add.scale.bar(x=0, y=-5)
          
        }
        
      }
      dev.off()
    }
    
  }
  
}


