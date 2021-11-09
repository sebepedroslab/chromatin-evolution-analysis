# functions and libraries
library(stringr)
library(readxl)

his_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")

# dictionary of modifications
dict_mods = data.frame(
    Rawmod =       c("Phospho", "Propy","3xMethyl","Trime","2xMethyl","Dime","Crotonyl","Acetyl"),
    Modification = c("Pho",     "Me1",  "Me3",     "Me3",  "Me2",     "Me2", "Cro",     "Ace")
)


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
mod_fn_list = list.files("../data-paired-searches/", pattern = "*.raw.csv", full.names = TRUE)
#mod_fn_list = mod_fn_list [ grepl("Gefoke", mod_fn_list )]

for (mod_fn in mod_fn_list) {
    
    # check if file is empty    
    mod_ref_info = file.info(mod_fn)
    if (mod_ref_info$size != 0 ) {
        
        # read ref seqs
        print(sprintf("Working on %s...", mod_fn))
        mod_ref = read.table(mod_fn, sep = ",", header = TRUE)
        
        if (grepl("Homo",mod_fn) & any(!grepl("^Hsap",mod_ref$MasterProteinAccessions)) ) {
            mod_ref$MasterProteinAccessions = paste("Hsap",mod_ref$MasterProteinAccessions, sep = "_")
        }
        
        # input files
        out_fn = str_remove(mod_fn, "\\.csv$")
        
        # drop duplicates and accessory names        
        mod_ref = mod_ref[!duplicated(mod_ref),]
        mod_ref$MasterProteinAccessions = gsub("; .*","",mod_ref$MasterProteinAccessions)
        
        # find previously seen histones
        mod_ref = merge(mod_ref, dic, all.x = T, all.y = F, by.x = "MasterProteinAccessions", by.y = "Protein_ID")
        
        # which histones are already annotated?
        is_annotated = grepl("[\\||_]Histone", mod_ref$MasterProteinAccessions)
        
        mod_ref$MasterProteinAccessions [ !is_annotated ] = paste(
            mod_ref$MasterProteinAccessions [ !is_annotated ], 
            "|Histone_",
            mod_ref$Histone_Type [ !is_annotated ],
            sep="")
        
        # keep histones
        mod_ref = mod_ref[grepl("[\\||_]Histone", mod_ref$MasterProteinAccessions), ]
        mod_ref = mod_ref[!grepl("[\\||_]Histone_NA", mod_ref$MasterProteinAccessions), ]
        
        if (nrow(mod_ref) > 0) {
            
            # homogenise histone annotations...
            mod_ref$MasterProteinAccessions = gsub("_Histone","|Histone",mod_ref$MasterProteinAccessions)
            mod_ref$MasterProteinAccessions = gsub("\\|HistoneH","|Histone_",mod_ref$MasterProteinAccessions)
            # discard unmodified peptides (pTMRS string)
            mod_ref = mod_ref[!is.na(mod_ref$ptmRSBestSiteProbabilities), ]
            mod_ref = mod_ref[mod_ref$ptmRSBestSiteProbabilities != "Too many isoforms", ]
            
            mod_sum = data.frame(
                Protein_ID = mod_ref$MasterProteinAccessions,
                Peptide_Sequence = mod_ref$Sequence,
                mods_string = mod_ref$ptmRSBestSiteProbabilities
            )
            
            # remove supplementary master proteins
            mod_sum$Protein_ID = str_split(mod_sum$Protein_ID, pattern = ";", simplify = T)[,1]
            
            mod_sua = data.frame()
            for (i in 1:nrow(mod_sum)) {
                
                pid = as.character(mod_sum$Protein_ID[i])
                his = str_remove(str_split(pid, pattern = "\\|")[[1]][2], "Histone_")
                pep = as.character(mod_sum$Peptide_Sequence[i])
                pep_pos = "res"
                mod_list = str_split(as.character(mod_sum$mods_string[i]), ";")[[1]]
                
                for (mod in mod_list) {
                    
                    mod = str_remove(mod, " ")
                    res = str_remove(str_split(mod, "\\(|\\)")[[1]][1] , pattern = "\\d+")
                    pos = as.numeric(str_remove(str_split(mod, "\\(|\\)")[[1]][1] , pattern = "[A-Za-z]*"))
                    typ = str_split(mod, "\\(|\\)")[[1]][2]
                    pro = as.numeric(str_split(mod, ":")[[1]][2])
                    
                    # print(paste(pid, his, pep, pep_pos, res, pos, typ, NA, pro))
                    mod_sui = data.frame(pid, his, pep, pep_pos, res, pos, typ, NA, pro)
                    mod_sua = rbind(mod_sua, mod_sui)
                    
                }
                
            }
            
            # header
            colnames(mod_sua) = c(
                "Protein_ID", "Histone_Type", "Peptide_Sequence",
                "Peptide_position", "Residue", "Position_in_peptide", "Rawmod",
                "Positions_in_protein", "ptmRS_Best_Site_Probabilities")
            
            # drop empty mods (Propionyl)
            mod_sua = mod_sua[mod_sua$Rawmod != "Propionyl",]
            # drop non-annotated histones
            mod_sua = mod_sua [ !is.na(mod_sua$Histone_Type) , ]
            
            # clean mod names
            mod_sua = merge(mod_sua, dict_mods, all.x=T, all.y=F, by.x="Rawmod", by.y="Rawmod")
            mod_sua = mod_sua[,c(
                "Protein_ID", "Histone_Type", "Peptide_Sequence",
                "Peptide_position", "Residue", "Position_in_peptide", "Modification",
                "Positions_in_protein", "ptmRS_Best_Site_Probabilities")
            ]
            
            # drop acetylations and empty mods
            mod_sua = mod_sua [ mod_sua$Modification != "Ace", ]
            mod_sua = mod_sua [ !is.na(mod_sua$Modification), ]
            
            if (nrow(mod_sua)>0) {
                
                
                # clean peptide
                pep_clean = mod_sua$Peptide_Sequence
                pep_clean = str_remove_all(pep_clean, pattern = "^\\[.*\\]\\.")
                pep_clean = str_remove_all(pep_clean, pattern = "\\.\\[.*\\]$")
                pep_clean = sub("\\s+$", "", gsub('(.{5})', '\\1 ', pep_clean))
                mod_sua$Clean_peptide = pep_clean
                
                # order by histone, prot, and peptide position
                start_peptide_pos = as.numeric(str_remove(str_split(mod_sua$Peptide_position, pattern = "-", simplify = T)[,1], "\\["))
                mod_sua = mod_sua[order(mod_sua$Histone_Type, mod_sua$Protein_ID,start_peptide_pos ,mod_sua$Position_in_peptide, mod_sua$Modification, -rank(mod_sua$ptmRS_Best_Site_Probabilities)),]
                
                # drop identical reports with low probability
                mod_sua_uniqid = unlist(lapply(1:nrow(mod_sua), function(i) {
                    paste(mod_sua[i,"Protein_ID"],mod_sua[i,"Clean_peptide"],mod_sua[i,"Modification"],mod_sua[i,"Position_in_peptide"], collapse = " ")
                }))
                mod_sux = mod_sua[ !duplicated( mod_sua_uniqid ) , ]
                
                # write.table(mod_sua, file = sprintf("%s.csv", out_fn), quote = F, sep = "\t", row.names = F)
                
                for (his in his_list) {
                    
                    mod_sui = mod_sua[mod_sua$Histone_Type == his,]
                    if (nrow(mod_sui)>0) {
                        mod_sui_f = data.frame(
                            Protein_ID = mod_sui$Protein_ID,
                            New_Protein_ID = mod_sui$Protein_ID,
                            Histone_Type = his,
                            is_variant = FALSE,
                            Clean_peptide = mod_sui$Clean_peptide,
                            Modification = mod_sui$Modification,
                            Residue = mod_sui$Residue,
                            Position_in_peptide = mod_sui$Position_in_peptide,
                            aa_in_reference = "",
                            Positions_in_Reference = "",
                            Positions_in_Protein = "",
                            ptmRS_Best_Site_Probabilities =mod_sui$ptmRS_Best_Site_Probabilities,
                            Peptide_Sequence = mod_sui$Peptide_Sequence,
                            Peptide_position =mod_sui$Peptide_position
                        )
                        
                        mod_sui_f = mod_sui_f[!duplicated(mod_sui_f),]
                        write.table(mod_sui_f, file = sprintf("%s-%s.csv", out_fn,his), quote = F, sep = "\t", row.names = F)
                    }
                    
                }
                
            } else {
                
                print(sprintf("Nothing to report for %s, file has no relevant modifications!", mod_fn))

            }
            
        } else {
            
            print(sprintf("Nothing to report for %s, file has modifications but not in histones!", mod_fn))
        }
        
    } else {
        
        print(sprintf("Nothing to report for %s, file is empty!", mod_fn))
        
    }
    
}



# system("for h in H2A macroH2A H2AZ H3 H4 H2B ; do for i in ../data-paired-searches/*-${h}.csv ;  do awk 'NR>1' $i ;  done > ../data-paired-searches/all_${h}.tsv ; done" )
