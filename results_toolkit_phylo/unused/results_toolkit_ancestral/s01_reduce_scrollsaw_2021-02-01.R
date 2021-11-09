#### Input ####

# load libraries
library(stringr)
library(seqinr)

clas_fn = "../data/gene_families_hmm.csv"
seqs_fn = "../results_toolkit_phylo/gene_sequences/"
outp_fn = "reduced_datasets"
tax_fn  = "../data/euk_taxonomy_annotated_2020-08-11.csv"


#### Define input ####

# load data
cla = read.table(clas_fn, header = F, sep = "\t", col.names = c("gene_class", "gene_type", "gene_fam", "domains", "search", "inflation", "min_size"))
tax = read.table(tax_fn, sep = "\t", header = T, stringsAsFactors = F)

# ignore complexes?
# cla = cla [ cla$gene_class != "Complexes", ]

# get list of gene families
gene_fams = as.character(cla$gene_fam)


# taxonomic groupings
scroll_sets = list(
  list(
    id = "euk",
    focus = "Macrogroup",
    tax = c( "Obazoa","Amoebozoa","CRuMs","Ancyromonadida","Malawimonadidae","ArchaeCry","SARHap","Hemimastigophora","Meteora","Discoba","Metamonada")
  ),
  list(
    id = "met",
    focus = "Group",
    tax = c( "Bilateria","Cnidaria","Placozoa","Ctenophora","Porifera")
  ),
  list(
    id = "opi",
    focus = "Group",
    tax = c( "Bilateria","Cnidaria","Placozoa","Ctenophora","Porifera","Choanoflagellata","Filasterea","Ichthyosporea","Corallochytrea","Fungi","Discicristoidea")
  ),
  list(
    id = "dia",
    focus = "Group",
    tax = c( "Streptophyta","Chlorophyta","Rhodophyta","Glaucophyta","Cryptista","Haptista","Stramenopiles","Alveolata","Rhizaria")
  )
)

scroll_sets = list(
  list(
    id = "eur",
    focus = "Macrogroup",
    tax = c( "Obazoa","Amoebozoa","CRuMs","Ancyromonadida","Malawimonadidae","ArchaeCry","SARHap","Hemimastigophora","Meteora","Discoba","Metamonada")
  )
)


#### Loop ####

# first through sets
# within each set, loop though gene families
# within each gene family, loop through 

for (set in scroll_sets) {
  
  # define set-specific vars
  set_foc = set$focus
  set_tax = set$tax
  set_id  = set$id
  
  for (fam in gene_fams) {
    
    # load sequences
    fas_list = list.files(pattern = sprintf("^euk.%s..*.fasta", fam), path = seqs_fn, full.names = T)
    fas = list()
    for (fai in fas_list) {
      fas_i = seqinr::read.fasta(fai)
      fas = c(fas, fas_i)
    }
    
    if (length(fas)>0){
      
      fas_sps = stringr::str_split(names(fas), pattern = "_", simplify = T)[,1]
      
      # empty list of genes to keep
      keep_genes = c()
      shared_genes = c()
      
      # start
      print(sprintf("%s | %s init genes %i", set_id, fam, length(fas_sps) ))
      # loop through set contents    
      for (ni in seq_along(set_tax)) {
        
        for (nj in seq_along(set_tax)) {
          
          if (ni < nj) {
            
            # define input      
            taxi = set_tax[ni]
            taxj = set_tax[nj]
            taxi_sps = tax[tax[,set_foc] == taxi, "Species" ]
            taxj_sps = tax[tax[,set_foc] == taxj, "Species" ]
            
            # define comparison-specific fasta files
            fasi = fas[fas_sps %in% taxi_sps]
            fasj = fas[fas_sps %in% taxj_sps]
            
            # proceed if there are sequences to align in both sets
            if (length(fasi) > 0 & length(fasj) > 0) {
              
              # write and align
              seqinr::write.fasta(fasi, file.out = sprintf("%s/tmp.i.fa", outp_fn), names = names(fasi))
              seqinr::write.fasta(fasj, file.out = sprintf("%s/tmp.j.fa", outp_fn), names = names(fasj))
              system(command = sprintf("diamond blastp --more-sensitive --max-target-seqs 2 -q %s/tmp.i.fa -d %s/tmp.j.fa -o %s/tmp.ij.csv --quiet --threads 6", outp_fn, outp_fn, outp_fn)  )
              system(command = sprintf("diamond blastp --more-sensitive --max-target-seqs 2 -q %s/tmp.j.fa -d %s/tmp.i.fa -o %s/tmp.ji.csv --quiet --threads 6", outp_fn, outp_fn, outp_fn)  )
              
              # load forward and backward
              # forward
              ali_ij = read.table(
                sprintf("%s/tmp.ij.csv", outp_fn, fam),
                col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
                sep="\t", header = F
              )
              # backward
              ali_ji = read.table(
                sprintf("%s/tmp.ji.csv", outp_fn, fam),
                col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"),
                sep="\t", header = F
              )
              
              # proceed if there are aligned sequences in both sets
              if (nrow(ali_ij)>0 & nrow(ali_ji)>0){
                
                # add species
                ali_ij$query_sps = stringr::str_split(ali_ij$qseqid, pattern = "_", simplify = T)[,1]
                ali_ij$subject_sps = stringr::str_split(ali_ij$sseqid, pattern = "_", simplify = T)[,1]
                ali_ji$query_sps = stringr::str_split(ali_ji$qseqid, pattern = "_", simplify = T)[,1]
                ali_ji$subject_sps = stringr::str_split(ali_ji$sseqid, pattern = "_", simplify = T)[,1]
                
                # add pair ids                
                ali_ij$pair = paste(ali_ij$qseqid, ali_ij$sseqid)
                ali_ji$pair = paste(ali_ji$sseqid, ali_ji$qseqid)
                
                # # keep only best hit for each gene pair
                # ali_ij = ali_ij[ !duplicated(ali_ij$pair), ]
                # ali_ji = ali_ji[ !duplicated(ali_ji$pair), ]
                
                # keep only best hit for each gene
                ali_ij = ali_ij[ !duplicated(ali_ij$qseqid), ]
                ali_ji = ali_ji[ !duplicated(ali_ji$qseqid), ]
                
                # remove poor hits
                ali_ij = ali_ij[ ali_ij$bitscore > 100, ]
                ali_ji = ali_ji[ ali_ji$bitscore > 100, ]
                
                # keep pairs present in both tables
                shared_pairs = ali_ij$pair[ali_ij$pair %in% ali_ji$pair]
                if (length(shared_pairs)>0){
                  shared_genes = c(
                    paste(taxi, ali_ij[ali_ij$pair %in% shared_pairs, "qseqid"], sep="|"),
                    paste(taxj, ali_ji[ali_ji$pair %in% shared_pairs, "qseqid"], sep="|")
                  )
                }
                
              } # end IF: check if there are alignments, skip if none
            } # end IF: check if there are sequences to align, skip if nothing
            
            # remove files
            system(command = sprintf("rm %s/tmp.*", outp_fn), ignore.stderr = T)
            
          } # end IF: which taxon sets to compare (to avoid redundancy)
          
          keep_genes = c(keep_genes,shared_genes)
          
        } # loop j along tax sets
      } # loop i along tax sets
      
      # report
      keep_genes = keep_genes[!is.na(keep_genes)]
      keep_genes_table = table(keep_genes)
      keep_genes_list  = names(keep_genes_table)[keep_genes_table > 1]
      print(sprintf("%s | %s kept genes %i", set_id, fam, length(keep_genes_list) ))
      write.table(
        data.frame(keep_genes_list), file=sprintf("%s/reduced.%s_%s.genes.txt", outp_fn, set_id, fam),
        row.names = F, quote = F, col.names = F)
      
      # save fasta 
      if (length(keep_genes_list)>0) { 
        fas_export_list = stringr::str_split(keep_genes_list, pattern = "\\|", simplify = T)[,2] 
      } else { 
        fas_export_list = keep_genes_list
      }
      write.fasta(
        fas[fas_export_list], file=sprintf("%s/reduced.%s_%s.genes.fasta", outp_fn, set_id, fam),
        names = keep_genes_list)
      
    } # end IF: are there any classified seqs at all?
  } # loop end: gene fam
} # loop end: scroll set

print("All done!")

