# functions and libraries
library(stringr)
library(stringi)
library(readxl)

his_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")
sps_list = read.table("../species_list_2021-02-10.txt", stringsAsFactors = FALSE)[,1]


# create dictionary of previously seen histone gene names
dic = data.frame()
for (his in his_list) { 
    
    print(sprintf("Find previously seen %s genes...", his))    
    dii = read_xlsx("../consensus_modifications_perseq.xlsx", sheet = his)
    dii = dii[,c("Protein_ID", "Histone_Type")]
    dii = dii[!duplicated(dii),]
    
    dic = rbind(dic, dii)
    
}

# clean
dic = dic [ !grepl(",", dic$Histone_Type), ]
dic = dic [ !grepl("_up_", dic$Protein_ID), ]
dic = dic [ !grepl("Uniprot", dic$Protein_ID), ]
dic = dic [ !is.na(dic$Histone_Type), ]
dic$Histone_Type = gsub(" .*","", dic$Histone_Type)
# drop homo species ids (not present in database)
# dic$Protein_ID = gsub("^Hsap_","",dic$Protein_ID)
# drop comments
dic$Protein_ID = gsub("\\|.*","",dic$Protein_ID)


#### Main loop ####

# list files to load
mod_fn_list = list.files("../data-reproducibility/", pattern = "*.replicates.csv", full.names = TRUE)
mod_fn_list = mod_fn_list [ !grepl("Methano", mod_fn_list) ]
mod_tot = data.frame()


for (mod_fn in mod_fn_list) {
    
    # check if file is empty    
    mod_ref_info = file.info(mod_fn)
    if (mod_ref_info$size != 0) {
        
        # read ref seqs
        print(sprintf("Working on %s...", mod_fn))
        mod_ref = read.table(mod_fn, sep = ",", colClasses = "character" , header = TRUE)
        colnames(mod_ref) = c("gene", "hptm_string", "hex_samples_clean",  "hex_samples", "PsmCount" , "Confidence", "sequence")
        
        # drop duplicates
        # mod_ref = mod_ref[!duplicated(mod_ref),]
        
        # keep histones
        if ( grepl("Homo",mod_fn) ) {
            mod_ref$gene [ !grepl("Hsap_", mod_ref$gene) ] = paste("Hsap", mod_ref$gene, sep = "_") [ !grepl("Hsap_", mod_ref$gene) ]
        } else if ( grepl("Tetrahy",mod_fn) )  {
            mod_ref$gene [ !grepl("Tthe_", mod_ref$gene) ] = paste("Tthe", mod_ref$gene, sep = "_") [ !grepl("Tthe_", mod_ref$gene) ]
        } else if ( grepl("Ectocarpus",mod_fn) )  {
            mod_ref$gene [ !grepl("Esil_", mod_ref$gene) ] = paste("Esil", mod_ref$gene, sep = "_") [ !grepl("Esil_", mod_ref$gene) ]
        } else if ( grepl("Arabidopsis",mod_fn) )  {
            mod_ref$gene [ !grepl("Atha_", mod_ref$gene) ] = paste("Atha", mod_ref$gene, sep = "_") [ !grepl("Atha_", mod_ref$gene) ]
        } else if ( grepl("dactilyum",mod_fn) )  {
            mod_ref$gene [ !grepl("Phatri_", mod_ref$gene) ] = paste("Phatri", mod_ref$gene, sep = "_") [ !grepl("Phatri_", mod_ref$gene) ]
        }
        mod_ref$gene_id = gsub("\\|.*","", mod_ref$gene)
        mod_ref$gene_id = gsub(";.*","", mod_ref$gene_id)

        mod_ref = mod_ref[ grepl("[\\||_]Histone", mod_ref[,1]) | mod_ref$gene_id %in% dic$Protein_ID, ]
        mod_ref = mod_ref[ !grepl("[\\||_]Histone_NA", mod_ref[,1]), ]

        
        # hex codes for peptide presence:  
        # 0400000001 green
        # 0000000001 gray (not found)
        # 0100000001 (peak found, blue)
        # 0300000001 (yellow, medium)
        
        # count number of samples
        mod_ref$samples_seen_list = lapply(mod_ref$hex_samples, function(v) as.numeric(grepl("^04|^03",  stringi::stri_sub(v, seq(1, stri_length(v), by = 10), length = 10))))
        
        mod_ref$samples_seen = unlist(lapply(mod_ref$hex_samples, function(v) sum(stringr::str_count(stringi::stri_sub(v, seq(1, stri_length(v), by = 10), length = 10), "^04|^03"))))
        # mod_ref$samples_seen = stringr::str_count(mod_ref$hex_samples, "04|03") #+ stringr::str_count(mod_ref$hex_samples, "2")
        mod_ref$samples_total = stringr::str_length(mod_ref$hex_samples)  / 10
        
        # concatenate
        run_id = gsub(".replicates.csv","",basename(mod_fn))
        mod_ref$run_id = run_id
        mod_tot = rbind(mod_tot, mod_ref)
        
    }
}

# add species
mod_tot$species = stringr::str_split(mod_tot$gene, "_", simplify = TRUE)[,1]
mod_tot$species = factor(mod_tot$species, levels = c(sps_list, "Hsap, replicate"))
mod_tot = mod_tot[!is.na(mod_tot$species), ]

# set apart the replicate double human experiment
mod_tot$species [ grepl("Homo-second", mod_tot$run_id) ] = "Hsap, replicate"

# keep interesting modifications
mod_tot_f = mod_tot[grepl("Acetyl\\b|Methyl\\b|Dimeyhyl\\b|Trimethyl\\b|Phospho\\b|Crotonyl\\b", mod_tot$hptm_string),]

# drop run id and remove duplicate peptides
# mod_tot$run_id = NULL
# mod_tot = mod_tot [!duplicated(mod_tot), ]



#### histograms ####

# plot per speceis
pdf("../replicability_hists.pdf", width = 10, height = 12)
layout(matrix(1:30, ncol=5))
for (spi in levels(mod_tot_f$species)) {
    
    mod_i = mod_tot_f [ mod_tot_f$species == spi, ]
    col = "lightblue2"
    if (spi %in% c("Phatri","Atha","Tthe","Esil","Hsap, replicate")) { col = "lightblue4" }
    if (nrow(mod_i) > 0) {

        mod_i_total = max(mod_i$samples_total)
        mod_i_per_single_peptide = unlist(lapply(1:nrow(mod_i), function(v) rep(mod_i$samples_seen[v], as.numeric(mod_i$PsmCount)[v])))
        # mod_i_per_single_peptide = mod_i$samples_seen
        
        # fraction seen in two samples
        seen_t = table(mod_i_per_single_peptide)
        rep_fra = sum(seen_t[names(seen_t) > 1]) / sum(seen_t)
        
        # plot
        hist(
            mod_i_per_single_peptide, 
            breaks = 0:max(c(mod_i_total, 12)), 
            col = col, 
            main = sprintf("%s\n%i modified PSMs (%.1fp in >1 sample)\n%i samples", spi, length(mod_i_per_single_peptide), rep_fra * 100, max(mod_i_total)),
            xlab = "", ylab = "",
            cex.axis = 0.8, cex.main = 0.8, las = 1
        )
        abline(v=max(mod_i_total), lty = 2, col = "red")
        
    } else {
        plot(0,0, col = NULL, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        text(0,0, sprintf("%s\nno data", spi), cex = 0.7)
    }
    
    
}


### Pie plots ###

# all
spi = "all hPTMs"
mod_i = mod_tot_f [ ! mod_tot_f$species %in% c("Phatri","Atha","Tthe","Esil","Hsap, replicate"),  ]

# fraction seen in various samples
mod_i_total = max(mod_i$samples_total)
mod_i_per_single_peptide = unlist(lapply(1:nrow(mod_i), function(v) rep(mod_i$samples_seen[v], as.numeric(mod_i$PsmCount)[v])))
# mod_i_per_single_peptide = mod_i$samples_seen
seen_t = table(mod_i_per_single_peptide)
seen_ti = seen_t [ 1:5 ]
seen_ex = sum(seen_t [ 6:length(seen_t) ])
if (seen_ex > 0) {
    seen_ti[">=6"] = seen_ex
}
pie(seen_ti, col = c("gray95","palevioletred1","palevioletred3","indianred1","indianred3","indianred4"), 
    labels = sprintf("%s | n=%i (%.1fp)",names(seen_ti), seen_ti, 100*seen_ti/sum(seen_ti)),
    main = sprintf("%s | n=%i\n(%.1fp in >1 sample)",spi, sum(seen_ti), 100*sum(seen_ti[2:length(seen_ti)]) / sum(seen_ti) ), cex.main = 0.8, cex= 0.6)



# all peptides with acetylations
spi = "all Ace"
mod_i = mod_tot_f [ ! mod_tot_f$species %in% c("Phatri","Atha","Tthe","Esil","Hsap, replicate"),  ]
mod_i = mod_i [ grepl("Acetyl",mod_i$hptm_string), ]

# fraction seen in various samples
mod_i_total = max(mod_i$samples_total)
mod_i_per_single_peptide = unlist(lapply(1:nrow(mod_i), function(v) rep(mod_i$samples_seen[v], as.numeric(mod_i$PsmCount)[v])))
# mod_i_per_single_peptide = mod_i$samples_seen
seen_t = table(mod_i_per_single_peptide)
seen_ti = seen_t [ 1:5 ]
seen_ex = sum(seen_t [ 6:length(seen_t) ])
if (seen_ex > 0) {
    seen_ti[">=6"] = seen_ex
}
pie(seen_ti, col = c("gray95","palevioletred1","palevioletred3","indianred1","indianred3","indianred4"), 
    labels = sprintf("%s | n=%i (%.1fp)",names(seen_ti), seen_ti, 100*seen_ti/sum(seen_ti)),
    main = sprintf("%s | n=%i\n(%.1fp in >1 sample)",spi, sum(seen_ti), 100*sum(seen_ti[2:length(seen_ti)]) / sum(seen_ti) ), cex.main = 0.8, cex= 0.6)



# all peptides with methylations
spi = "all Me1/2/3"
mod_i = mod_tot_f [ ! mod_tot_f$species %in% c("Phatri","Atha","Tthe","Esil","Hsap, replicate"),  ]
mod_i = mod_i [ grepl("Methyl|Dimethyl|Trimethyl",mod_i$hptm_string), ]

# fraction seen in various samples
mod_i_total = max(mod_i$samples_total)
mod_i_per_single_peptide = unlist(lapply(1:nrow(mod_i), function(v) rep(mod_i$samples_seen[v], as.numeric(mod_i$PsmCount)[v])))
# mod_i_per_single_peptide = mod_i$samples_seen
seen_t = table(mod_i_per_single_peptide)
seen_ti = seen_t [ 1:5 ]
seen_ex = sum(seen_t [ 6:length(seen_t) ])
if (seen_ex > 0) {
    seen_ti[">=6"] = seen_ex
}
pie(seen_ti, col = c("gray95","palevioletred1","palevioletred3","indianred1","indianred3","indianred4"), 
    labels = sprintf("%s | n=%i (%.1fp)",names(seen_ti), seen_ti, 100*seen_ti/sum(seen_ti)),
    main = sprintf("%s | n=%i\n(%.1fp in >1 sample)",spi, sum(seen_ti), 100*sum(seen_ti[2:length(seen_ti)]) / sum(seen_ti) ), cex.main = 0.8, cex= 0.6)

dev.off()




# how many peptides are shared between replicates?
# spi = "Hsap"
# mod_i = mod_tot_f [ mod_tot_f$species == "Hsap",  ]
# mod_j = mod_tot_f [ mod_tot_f$species == "Hsap, replicate",  ]
# 
# sum(unique(mod_i$sequence) %in% unique(mod_j$sequence)) / length(unique(mod_i$sequence))
# sum(unique(mod_j$sequence) %in% unique(mod_i$sequence)) / length(unique(mod_j$sequence))
# 
