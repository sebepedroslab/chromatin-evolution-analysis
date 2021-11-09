#### Input ####

# load libraries
library(ape)
library(stringr)

clas_fn = "../data/gene_families_hmm.csv"
alig_fn = "alignments/"
outp_fn = "trees_reduced"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"


#### Define input ####

# load data
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

# ignore complexes
cla = cla [ cla$gene_class != "Complexes", ]

# get list of gene families
gene_fams = as.character(cla$gene_fam)


# taxonomic groupings
scroll_sets = list(
  list(
    id = "euk",
    focus = "Macrogroup", 
    tax = c( "Obazoa","Amoebozoa","CRuMs","Ancyromonadida","Malawimonadidae","ArchaeCry","SARHap","Hemimastigophora","Meteora","Discoba","Metamonada")
  )
  # list(
  #   id = "met",
  #   focus = "Group", 
  #   tax = c( "Bilateria","Cnidaria","Placozoa","Ctenophora","Porifera")
  # ),
  # list(
  #   id = "opi",
  #   focus = "Group", 
  #   tax = c( "Bilateria","Cnidaria","Placozoa","Ctenophora","Porifera","Choanoflagellata","Filasterea","Ichthyosporea","Corallochytrea","Fungi","Discicristoidea")
  # ),
  # list(
  #   id = "pla",
  #   focus = "Group", 
  #   tax = c( "Streptophyta","Chlorophyta","Rhodophyta","Glaucophyta","Cryptista","Haptista","Stramenopiles","Alveolata","Rhizaria")
  # )
)



#### Loop ####

# first through sets
# within each set, loop though gene families
# within each gene family, loop through 

for (fam in gene_fams) {

  # load alignments
  ali_fn = sprintf("",alig_fn, fam)
  ali = read.table(
    sprintf("%s/euk.%s.diamond.csv.gz", alig_fn, fam),
    col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
    sep="\t", header = F
  )
  ali$query_sps = stringr::str_split(ali$qseqid, pattern = "_", simplify = T)[,1]
  ali$subject_sps = stringr::str_split(ali$sseqid, pattern = "_", simplify = T)[,1]
  
  for (set in scroll_sets) {
    
    # define set-specific vars
    set_foc = set$focus
    set_tax = set$tax
    set_id  = set$id
    
    # empty list of genes to keep
    keep_genes = c()
	shared_genes = c()
    
    # loop through set contents    
    for (ni in seq_along(set_tax)) {
      
      for (nj in seq_along(set_tax)) {
        
        if (ni < nj) {
          
          # define input      
          taxi = set_tax[ni]
          taxj = set_tax[nj]
          taxi_sps = tax[tax[,set_foc] == taxi, "Species" ]
          taxj_sps = tax[tax[,set_foc] == taxj, "Species" ]
          
          # subset forward and backward
          ali_ij = ali [ ali$query_sps %in% taxi_sps & ali$subject_sps %in% taxj_sps ,]
          ali_ji = ali [ ali$query_sps %in% taxj_sps & ali$subject_sps %in% taxi_sps ,]
          
          # # keep only best hit
          ali_ij = ali_ij[ !duplicated(ali_ij$qseqid), ]
          ali_ji = ali_ji[ !duplicated(ali_ji$qseqid), ]
          
          # keep pairs present in both tables
          ali_ij$pair = paste(ali_ij$qseqid, ali_ij$sseqid)
          ali_ji$pair = paste(ali_ji$sseqid, ali_ji$qseqid)
          shared_pairs = ali_ij$pair[ali_ij$pair %in% ali_ji$pair]
          shared_genes = c(
            ali_ij[ali_ij$pair %in% shared_pairs, "qseqid"],
            ali_ji[ali_ji$pair %in% shared_pairs, "qseqid"]
          )
          
        }

		# store genes
	    keep_genes = unique(c(keep_genes,shared_genes))

        
      }
      
    }
    
    # report  
    print(sprintf("%s | %s init genes %i", set_id, fam, length(unique(c(ali$sseqid, ali$qseqid)))))
    print(sprintf("%s | %s kept genes %i", set_id, fam, length(keep_genes)))
    
    write.table(
      data.frame(keep_genes), file=sprintf("%s/reduced.%s_%s.genes.txt", alig_fn, set_id, fam),
      row.names = F, quote = F, col.names = F)
    
  }
  
}

print("All done!")

