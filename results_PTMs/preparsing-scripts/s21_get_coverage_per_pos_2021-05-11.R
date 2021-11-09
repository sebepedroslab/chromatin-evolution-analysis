# functions and libraries
library(stringr)
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
dic = dic [ !grepl("_up_", dic$Histone_Type), ]
dic = dic [ !is.na(dic$Histone_Type), ]
# drop homo species ids (not present in database)
# dic$Protein_ID = gsub("^Hsap_","",dic$Protein_ID)
# drop comments
dic$Protein_ID = gsub("\\|.*","",dic$Protein_ID)

#### Main loop ####

# list files to load
mod_fn_list = list.files("../data-coverage/", pattern = "*.perpos_cov_seqids.csv", full.names = TRUE)
mod_tot = data.frame()

for (mod_fn in mod_fn_list) {
    
    cop_fn = sub(".perpos_cov_seqids.csv", ".perpos_cov_peptides.csv", mod_fn)
    
    # check if file is empty    
    mod_ref_info = file.info(mod_fn)
    cop_ref_info = file.info(cop_fn)
    if (mod_ref_info$size != 0 & cop_ref_info$size != 0) {
        
        # read ref seqs
        print(sprintf("Working on %s...", mod_fn))
        mod_ref = read.table(mod_fn, sep = ",", header = TRUE)
        cop_ref = read.table(cop_fn, sep = ",", header = TRUE, col.names = c("TargetProteinsUniqueSequenceID", "SequencesinProtein", "PositionsinProtein"))

        if (grepl("Homo",mod_fn) & any(!grepl("^Hsap",mod_ref$Accession)) ) {
            mod_ref$Accession = paste("Hsap",mod_ref$Accession, sep = "_")
        }
        
        # drop duplicates and accessory names        
        mod_ref = mod_ref[!duplicated(mod_ref),]
        
        # find previously seen histones
        mod_ref = merge(mod_ref, dic, all.x = TRUE, all.y = FALSE, by.x = "Accession", by.y = "Protein_ID")
        
        # which histones are already annotated?
        is_annotated = grepl("[\\||_]Histone", mod_ref$Accession)
        
        mod_ref$Accession [ !is_annotated ] = paste(
            mod_ref$Accession [ !is_annotated ], 
            "|Histone_",
            mod_ref$Histone_Type [ !is_annotated ],
            sep="")
        
        # keep histones
        mod_ref = mod_ref[grepl("[\\||_]Histone", mod_ref$Accession), ]
        mod_ref = mod_ref[!grepl("[\\||_]Histone_NA", mod_ref$Accession), ]
        
        
        cop_ref = cop_ref [ cop_ref$PositionsinProtein != "", ]
        cop_ref = cop_ref[!duplicated(cop_ref),]
        cop_ref$start = as.numeric(stringr::str_split(gsub("\\[|\\]","",cop_ref$PositionsinProtein), "-", simplify = TRUE)[,1])
        cop_ref$stop = as.numeric(stringr::str_split(gsub("\\[|\\]","",cop_ref$PositionsinProtein), "-", simplify = TRUE)[,2])
        
        cop_mer = merge(mod_ref, cop_ref, by.x = "UniqueSequenceID", by.y = "TargetProteinsUniqueSequenceID", all.x = TRUE, all.y = FALSE)
        
        
        
        # store
        mod_tot = rbind(mod_tot, cop_mer)
        
    }
    
}

# sort by protein and coverage
mod_tot = mod_tot [ order(mod_tot$Accession, mod_tot$start, decreasing = FALSE),  ]

# keep highest coverage estimate per prot
mod_tot = mod_tot [ !is.na(mod_tot$Accession), ]
mod_tot$species = stringr::str_split(mod_tot$Accession, "_", simplify = TRUE)[,1]
mod_tot$species = factor(mod_tot$species, levels = sps_list)
mod_tot = mod_tot[ !is.na(mod_tot$species), ]
mod_tot = mod_tot [ order(mod_tot$species), ]
mod_tot$SequencesinProtein  = gsub("\\.", "", mod_tot$SequencesinProtein)

# get histone type
mod_tot$Histone_Type [ is.na(mod_tot$Histone_Type) ] = sub("Histone_","",stringr::str_split(mod_tot$Accession [ is.na(mod_tot$Histone_Type) ], "\\|", simplify = TRUE)[,2])

write.table(mod_tot, file="../coverage_per_position.csv", quote = FALSE, sep = "\t", row.names = FALSE)
