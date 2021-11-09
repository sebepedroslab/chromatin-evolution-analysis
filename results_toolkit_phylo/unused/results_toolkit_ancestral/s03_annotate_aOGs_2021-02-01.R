#### Input ####

# load libraries
library(stringr)
library(alluvial)
library(ape)

clas_fn = "../data/gene_families_hmm.csv"
aogs_fn = "reduced_trees/"
dat_fn  = "../results_toolkit_phylo/orthogroups_euk.csv"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"
outp_fn = "."
phyl_fn = "../results_toolkit_phylo/data/species_tree.newick"


#### Define input ####

# load data
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)
dat = read.table(dat_fn, header = T, stringsAsFactors = F)
dat$gene_class = stringr::str_split(dat$orthogroup, pattern = "\\.", simplify = T)[,1]

# read species tree
phyl = ape::read.tree(file = phyl_fn)
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



# which ancestral set to analyse?
scroll_sets = list(
  list(set = "euk", set_tax = "Eukaryota"),
  list(set = "met", set_tax = "Metazoa"),
  list(set = "opi", set_tax = "Opisthokonta"),
  list(set = "pla", set_tax = "ArchaeCry")
)
scroll_sets = list(
  list(set = "eur", set_tax = "Eukaryota")
)

#### Loop scrollsaw sets ####

# get list of gene families
gene_fams = as.character(cla$gene_fam)

for (seti in scroll_sets){ 
  
  # get set-specific info
  set = seti$set
  set_tax = seti$set_tax
  
  # create set-specific dataframe
  aog_all = data.frame()
  
  # restrict presence table to taxa that directly descend from ancestor of interest
  if (set_tax == "Eukaryota") {
    root_descendants = as.character(phyl_edge[phyl_edge$edge_start == length(phyl$tip.label)+1,"taxa"])
  } else {
    root_descendants = as.character(phyl_edge[phyl_edge$edge_start == phyl_edge[phyl_edge$taxa == set_tax,"edge_end"],"taxa"])
  }
  root_descendants = root_descendants[!is.na(root_descendants)]
  root_descendants_list = list()
  phyl_subtrees = subtrees(phyl)
  # get species (tips) descending from each LECA descendant subtree
  for (dei in root_descendants) {
    if (dei %in% phyl$node.label) {
      dei_subtree = phyl_subtrees[[ which(phyl$node.label == dei) ]]
      root_descendants_list[[dei]] = dei_subtree$tip.label
    } else if (dei %in% phyl$tip.label) {
      root_descendants_list[[dei]] = dei
    } else {
      root_descendants_list[[dei]] = NA
    }
  }
  root_descendants_table = data.frame(
    taxon_group = stringr::str_remove(names(unlist(root_descendants_list)), "\\d+"), 
    species = unlist(root_descendants_list))
  
  #### Loop gene families ####
  # gene_fams = "Hist_deacetyl"  
  # for each gene family:
  pdf(sprintf("%s/ancestral_OGs.%s.class.pdf", outp_fn, set),height = 8, width = 6)
  for (fam in gene_fams) {
    
    print(sprintf("%s | %s",set,fam))
    if (file.exists(sprintf("%s/reduced.%s_%s.possom.ortholog_groups.csv", aogs_fn, set, fam))) {
      
      # read ancestral orthologs
      aog = read.table(sprintf("%s/reduced.%s_%s.possom.ortholog_groups.csv", aogs_fn, set, fam), sep="\t", header = T)
      colnames(aog) = c("gene_wtax","aOG","aOG_to")
      aog$gene = stringr::str_split(aog$gene_wtax, pattern="\\|", simplify = T)[,2]
      aog$species = stringr::str_split(aog$gene, pattern="_", simplify = T)[,1]
      # get reference orthologs
      ogi = dat[dat$gene_class == fam,]
      aog = merge(aog, ogi, by.x = "gene", by.y = "gene", all.x = T, all.y = F)
      
      # get ref OG label
      is_like_og = grepl(":like:",aog$orthogroup)
      aog$orthogroup_label = NA
      if (any(is_like_og) & any(!is_like_og)) {
        aog[is_like_og,"orthogroup_label"] =  paste("like",stringr::str_split(aog[is_like_og,"orthogroup"],  pattern = ":", simplify = T)[,3], sep=":")
        aog[!is_like_og,"orthogroup_label"] = stringr::str_split(aog[!is_like_og,"orthogroup"], pattern = ":", simplify = T)[,2]
      } else if (any(is_like_og) & !any(!is_like_og)) {
        aog[is_like_og,"orthogroup_label"] =  paste("like",stringr::str_split(aog[is_like_og,"orthogroup"],  pattern = ":", simplify = T)[,3], sep=":")
      } else if (!any(is_like_og) & any(!is_like_og)) {
        aog[,"orthogroup_label"] = stringr::str_split(aog[,"orthogroup"], pattern = ":", simplify = T)[,2]
      }
      
      # stratify OG according to strength of evidence in the ancestor
      aog = merge(aog, root_descendants_table, by.x = "species", by.y = "species", all.x = T , all.y = F)
      aog_agg = aggregate(aog$taxon_group, by=list(aOG = aog$aOG), function(x) paste(unique(sort(x)), collapse = ","))
      aog = merge(aog, aog_agg, by.x = "aOG", by.y = "aOG", all.x = T, all.y = F)  
      
      # store ancestral OG info
      colnames(aog_agg) = c("aOG","ancestral_evidence_string")
      aog_agg$gene_fam = fam
      aog_agg$num_lineages = aggregate(aog$taxon_group, by=list(aOG = aog$aOG), function(x) length(unique(sort(x))))[,2]
      aog_agg$num_species = aggregate(aog$species, by=list(aOG = aog$aOG), function(x) length(unique(x)))[,2]
      aog_agg$num_genes = aggregate(aog$gene, by=list(aOG = aog$aOG), function(x) length(x))[,2]
      aog_agg$refOGs = aggregate(aog$orthogroup, by=list(aOG = aog$aOG), function(x) paste(unique(sort(x)), collapse = ","))[,2]
      aog_agg$refOGs_label = aggregate(aog$orthogroup_label, by=list(aOG = aog$aOG), function(x) paste(unique(sort(x)), collapse = ","))[,2]
      aog_agg$list_genes = aggregate(aog$gene, by=list(aOG = aog$aOG), function(x) paste(unique(sort(x)), collapse = ","))[,2]
      
      # discard non-ancestral aOGs
      aog = aog[grepl(aog$x,pattern = ","),]
      # aog_agg = aog_agg[aog_agg$num_lineages>1,]
      
      if (nrow(aog)>0){
        # visualise classification in alluvial plot
        aog_alluvial = as.data.frame(table(aog$aOG, aog$orthogroup))
        aog_alluvial = aog_alluvial[aog_alluvial$Freq > 2,]
        if (nrow(aog_alluvial)>1){
          colnames(aog_alluvial) = c("ancestral_OG", "OG", "Freq")
          alluvial(
            aog_alluvial[,1:2], 
            freq = aog_alluvial[,3], 
            col = c("lightblue"), gap.width = 0.1, blocks = T, 
            cex=0.4, border=c(NA))
          title(sprintf("ancestral OGs %s %s", set, fam),
                sub = sprintf("n genes = %i\nn aOG = %i", sum(aog_alluvial$Freq), length(unique(aog_alluvial$ancestral_OG))))
        } else {
          plot(0,0, xlab = "", ylab = "", frame.plot = F, axes = F, col="blue", pch=19)
          title(sprintf("ancestral OGs %s %s", set, fam),
                sub = paste(as.character(aog_alluvial[,1]),"=", as.character(aog_alluvial[,2]), "\nn genes =", aog_alluvial[,3]))
        }
      }
      
      # store data for this family
      aog_all = rbind(aog_all, aog_agg)
      
    }
  }
  dev.off()
  
  # save ancestor-specific dataframe
  write.table(aog_all, file=sprintf("%s/ancestral_OGs.%s.class.csv", outp_fn, set),
              quote = F, sep="\t", row.names = F)
  
}

