# functions and libraries
library(stringr)
library(readxl)

his_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")

# dictionary of modifications
dict_mods = data.frame(
    Rawmod =       c("Phospho", "Propy","3xMethyl","Trime","2xMethyl","Dime","Crotonyl","Acetyl"),
    Modification = c("Pho",     "Me1",  "Me3",     "Me3",  "Me2",     "Me2", "Cro",     "Ace")
)

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
mod_fn_list = list.files("../data-coverage/", pattern = "*.coverage.csv", full.names = TRUE)
mod_tot = data.frame()

for (mod_fn in mod_fn_list) {
    

    # check if file is empty    
    mod_ref_info = file.info(mod_fn)
    if (mod_ref_info$size != 0 ) {
        
        # read ref seqs
        print(sprintf("Working on %s...", mod_fn))
        mod_ref = read.table(mod_fn, sep = ",", header = TRUE)
        
        if (grepl("Homo",mod_fn) & any(!grepl("^Hsap",mod_ref$Accession)) ) {
            mod_ref$Accession = paste("Hsap",mod_ref$Accession, sep = "_")
        }
        
        # input files
        out_fn = str_remove(mod_fn, "\\.csv$")
        
        # drop duplicates and accessory names        
        mod_ref = mod_ref[!duplicated(mod_ref),]
        
        # find previously seen histones
        mod_ref$Accession = gsub("\\|.*","",mod_ref$Accession)
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
        
        # store
        mod_tot = rbind(mod_tot, mod_ref)
        
    }
    
}

# sort by protein and coverage
mod_tot = mod_tot [ order(mod_tot$Accession, mod_tot$Coverage, decreasing = TRUE),  ]

# keep highest coverage estimate per prot
mod_tot = mod_tot [ !duplicated(mod_tot$Accession), ]
mod_tot$species = stringr::str_split(mod_tot$Accession, "_", simplify = TRUE)[,1]
mod_tot$species = factor(mod_tot$species, levels = sps_list)
mod_tot = mod_tot[ !is.na(mod_tot$species), ]
mod_tot = mod_tot [ order(mod_tot$species), ]
# get histone type
# mod_tot$Histone_Type [ is.na(mod_tot$Histone_Type) ] = sub("Histone_","",stringr::str_split(mod_tot$Accession [ is.na(mod_tot$Histone_Type) ], "\\|", simplify = TRUE)[,2])
mod_tot$Histone_Type [ is.na(mod_tot$Histone_Type) ] = ""


write.table(mod_tot, file="../coverage_per_histone.csv", quote = FALSE, sep = "\t", row.names = FALSE)

pdf("../coverage_per_histone.max_value.pdf", width = 3.5, height = 5)
layout(matrix(1:3, ncol = 3))
for (his in his_list) {
    
    ### Coverage ###
    modi = mod_tot [ mod_tot$Histone_Type == his | grepl(his, mod_tot$Accession), ]
    moda = aggregate(data = modi, Coverage ~ species, max)
    modt = merge(moda, data.frame(species = sps_list), by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
    b=barplot(rev(modt$Coverage), names = rev(modt$species), horiz = TRUE, las = 1, xlim = c(0,100), main = his, cex.names =0.5, border = NA, xlab = "coverage %")
    text(x=0, b, labels = rev(sprintf("%.1f",modt$Coverage)),col="darkred",pos=4,cex=0.5)
    title(sub = sprintf("mean=%.2f,median=%.2f", mean(rev(modt$Coverage), na.rm = T), median(rev(modt$Coverage), na.rm = T)), cex.sub = 0.5)
    
    ### Num modifications ###
    mod = read_excel("../consensus_modifications_perseq.xlsx", sheet=his)
    mou = read_excel("../consensus_modifications_perseq.xlsx", sheet=paste("ub",his, sep = ""))
    mod = rbind(mod, mou)
    mod = mod[!is.na(mod$Positions_in_Reference),]
    # add species
    mod$species = stringr::str_split(mod$Protein_ID, pattern = "_", simplify = TRUE)[,1]
    mod$species = factor(mod$species, levels = sps_list)
    # unique modifications
    mod$Modification = factor(mod$Modification, levels = c("Ace","Me1","Me2","Me3","Cro","Pho","Ubi"))
    # unique modifications
    mof = mod [ , c("species", "Positions_in_Reference" , "Modification" )]
    mof = unique(mof)
    mot = table(mof$species)
    b=barplot(rev(mot), names = rev(names(mot)), horiz = TRUE, las = 1, xlim = c(0,40), main = his, cex.names =0.5, border = NA, xlab = "hPTMs")
    text(x=0, b, labels = rev(sprintf("%i",mot)),col="darkred",pos=4,cex=0.5)
    title(sub = sprintf("mean=%.2f,median=%.2f", mean(rev(mot), na.rm = T), median(rev(mot), na.rm = T)), cex.sub = 0.5)

    
    ### Modifications/covered position ###
    modw = aggregate(data = modi, Coverage ~ species, which.max)
    modw_ixs = modw$Coverage
    modl = aggregate(data = modi, NumberOfAminoAcids ~ species, c)
    modl_naa = modl$NumberOfAminoAcids
    modl_naa_maxcov = unlist(lapply(1:length(modl_naa), function(i) { 
        modl_naa[[i]] [ modw_ixs[i] ]
    }))
    
    modtm = merge(modt, data.frame(species = modl$species, num_positions = modl_naa_maxcov), by = "species", all.x = TRUE)
    modtm$covered_positions = modtm$Coverage / 100 * modtm$num_positions
    b=barplot(rev(mot) / rev(modtm$covered_positions), names = rev(names(mot)), horiz = TRUE, las = 1, main = his, cex.names =0.5, border = NA, xlab = "hPTMs/covered pos", xlim = c(0,1))
    text(x=0, b, labels = sprintf("%.2f", rev(mot) / rev(modtm$covered_positions)),col="darkred",pos=4,cex=0.5)
    title(sub = sprintf("mean=%.2f,median=%.2f", mean(rev(mot) / rev(modtm$covered_positions), na.rm = T), median(rev(mot) / rev(modtm$covered_positions), na.rm = T)), cex.sub = 0.5)
    
}
dev.off()
