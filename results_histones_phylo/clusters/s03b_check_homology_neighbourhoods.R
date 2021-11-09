# libraries
library(stringr)
library(ape)
library(dplyr)
library(zoo)

# input
gen_fo = "/home/xavi/Documents/Lab/prova-genomes/"
sps_fn = "../../data//validation_species_map.csv"
ali_fo = "/home/xavi/Documents/Lab/prova-genomes/green/OrthoFinder/Results_Mar14/WorkingDirectory/"
aix_fn = "/home/xavi/Documents/Lab/prova-genomes/green/OrthoFinder/Results_Mar14/WorkingDirectory/SpeciesIDs.txt"
gix_fn = "/home/xavi/Documents/Lab/prova-genomes/green/OrthoFinder/Results_Mar14/WorkingDirectory/SequenceIDs.txt"
pai_fn = "summary_pairs.csv"

# get species list
sps = read.table(sps_fn, header = T, stringsAsFactors = F)
sps_list = sps$species_assembly
sps_list = c("Crei","Vcar","Caulen","Cvar","Cocsub","Spun","Bden","Ttra","Klenit")
sps_list = c("Crei","Vcar","Cvar")

# get sequences and species indexes
aix = read.table(aix_fn, sep = ":", col.names = c("orthofinder_spsid", "file"))
aix$file = gsub(" ", "", aix$file)
aix$species = stringr::str_split(aix$file,"_", simplify = T)[,1]
gix = read.table(gix_fn, sep=":", col.names = c("orthofinder_geneid","gene"))
gix$gene = gsub(" ", "", gix$gene)

# get pairs of interest
pai = read.table(pai_fn, header = T, sep = "\t")
pai_h2ah2b = pai[ pai$pair_hists_ori == "H2A,H2B hh", ]
pai_h3h4 = pai[ pai$pair_hists_ori == "H3,H4 hh", ]

# get gene indexes from gtf; FUNCTION
gene_index_from_gtf = function(gtf_fn) {
    
    # load GTF
    gti = read.gff(gtf_fn)
    gti = gti[gti$type == "transcript",]
    # ensure order
    gti = gti[order(gti$seqid, gti$start), ]
    
    # get scaffold numeric id
    gei = data.frame(
        gene_original = gsub("transcript_id |\"", "", stringr::str_split(gti$attributes, ";", simplify = T)[,1]),
        chr_num = as.numeric(gti$seqid)
    )
    # add gene numbering
    gei = gei %>%
        dplyr::group_by(chr_num) %>%
        dplyr::mutate(gen_num = 1:n())
    gei$species = stringr::str_split(gei$gene_original, "_", simplify = T)[,1]
    
    return(gei)
    
}


#### Work: get pairs ####

for (spi in sps_list[1]) {
    
    # define input    
    spi_a = spi
    spi_h = sps[sps$species_assembly == spi, "species_histonome"]
    spi_ix = aix[aix$species == spi, "orthofinder_spsid"]
    print(sprintf("Get histone pairs from %s (syn=%s)", spi_a, spi_h))
    
    # get gene index
    spi_gtf_fn = sprintf("%s/%s_long.annot.gtf", gen_fo, spi_a)
    gei = gene_index_from_gtf(spi_gtf_fn)
    gei = merge(gei, gix, by.x = "gene_original", by.y = "gene", all.x = T, all.y = F)
    
    for (spj in sps_list) {
        
        if (spi != spj) {
            
            # define input    
            spj_a = spj
            spj_h = sps[sps$species_assembly == spj, "species_histonome"]
            spj_ix = aix[aix$species == spj, "orthofinder_spsid"]
            
            # get gene index            
            spj_gtf_fn = sprintf("%s/%s_long.annot.gtf", gen_fo, spj_a)
            gej = gene_index_from_gtf(spj_gtf_fn)
            gej = merge(gej, gix, by.x = "gene_original", by.y = "gene", all.x = T, all.y = F)
            
            # get pairwise alignments
            ali_ij = read.table(
                sprintf("%s/Blast%i_%i.txt.gz", ali_fo, spi_ix, spj_ix),
                col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))
            
            # rolling function along species i
            zoo::rollapply(gei, width=10, function(x) {
                
                
                
            })
            
            
        }
        
    }
    
}
