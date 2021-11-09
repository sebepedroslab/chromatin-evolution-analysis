# libraries
library(stringr)
library(ape)
library(dplyr)

# input
gen_fo = "/users/asebe/xgraubove/genomes/data/"
sps_fn = "../../data//validation_species_map.csv"

sps = read.table(sps_fn, header = T, stringsAsFactors = F)
#sps_list = c("Crei","Vcar","Caulen","Cvar","Cocsub","Spun","Bden","Ttra","Klenit")
sps_list = c("Tvag","Ehux","Honfer","Sappar","Alblai","Ngad","Phatri","Aano","Porpur","Porumb","Chocri","Gracho","Mpus","Otau","Cvar","Caulen","Vcar","Crei","Klenit","Sbic","Ttra","Rall","Cang","Amac","Spun","Bden","Pisp","Gpro","Coerev","Rirr","Morelo","Umay","Ccin","Cneo","Scer","Tmel","Ncra","Spom","Cowc","Emue","Tetwil","Morvir","Hvul","Spis","Adig","Exapal","Nvec","Emul","Phoaus","Hduj","Cscu","Lpol","Spur","Skow","Drer","Hsap")


for (spi in sps_list) {

    # define input    
    spa = spi
    sph = sps[sps$species_assembly == spi, "species_histonome"]
    print(sprintf("Prepare %s (syn=%s)", spa, sph))
    
    if (file.exists( sprintf("%s/%s_long.annot.gtf", gen_fo, spa) )) {
    # load GTF
    gtf = read.gff(sprintf("%s/%s_long.annot.gtf", gen_fo, spa))
    gtf = gtf[gtf$type == "transcript",]
    # ensure order
    gtf = gtf[order(gtf$seqid, gtf$start), ]
    
    # get scaffold numeric id
    gen = data.frame(
        gene_original = gsub("transcript_id |\"", "", stringr::str_split(gtf$attributes, ";", simplify = T)[,1]),
        chr_num = as.numeric(gtf$seqid)
    )
    # add gene numbering
    gen = gen %>%
        dplyr::group_by(chr_num) %>%
        dplyr::mutate(gen_num = 1:n())
    gen$species = stringr::str_split(gen$gene_original, "_", simplify = T)[,1]
    
    # compose final sequence name
    gen$gene_evolclust = paste(gen$species, gen$chr_num, gen$gen_num, sep="_")

    write.table(gen[,c("gene_original","gene_evolclust")], sprintf("data-evolclust/%s.gene_dict.csv", spa), sep = "\t", quote = F, row.names = F, col.names = F)
    } 
}
