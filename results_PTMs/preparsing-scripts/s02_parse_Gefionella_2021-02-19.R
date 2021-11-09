# functions and libraries
library(readxl)
library(stringr)

# input files
mod_fn = "../data/Gefoke-allfiles-allmod_PSMs.xlsx"
out_fn = str_remove(mod_fn, "\\.xlsx$")
his_list = c("H3","H4","H2A","H2B","macroH2A","H2AZ")

# dictionary of modifications
dict_mods = data.frame(
    Rawmod =       c("Phospho", "Propy","3xMethyl","2xMethyl","Crotonyl","Acetyl"),
    Modification = c("Pho",     "Me1",  "Me3",     "Me2",     "Cro",     "Ace")
)

# read ref seqs
mod_ref = read_excel(mod_fn, sheet="Full1")


# keep histones
mod_ref = mod_ref[grepl("\\|Histone", mod_ref$`Master Protein Accessions`), ]
# discard unmodified peptides (pTMRS string)
mod_ref = mod_ref[!is.na(mod_ref$`ptmRS: Best Site Probabilities`), ]
mod_ref = mod_ref[mod_ref$`ptmRS: Best Site Probabilities` != "Too many isoforms", ]

mod_sum = data.frame(
    Protein_ID = mod_ref$`Master Protein Accessions`,
    Peptide_Sequence = mod_ref$`Annotated Sequence`,
    mods_string = mod_ref$`ptmRS: Best Site Probabilities`
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
        
        print(paste(pid, his, pep, pep_pos, res, pos, typ, NA, pro))
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

# clean mod names
mod_sua = merge(mod_sua, dict_mods, all.x=T, all.y=F, by.x="Rawmod", by.y="Rawmod")
mod_sua = mod_sua[,c(
    "Protein_ID", "Histone_Type", "Peptide_Sequence",
    "Peptide_position", "Residue", "Position_in_peptide", "Modification",
    "Positions_in_protein", "ptmRS_Best_Site_Probabilities")
]

# clean peptide
pep_clean = mod_sua$Peptide_Sequence
pep_clean = str_remove_all(pep_clean, pattern = "^\\[.*\\]\\.")
pep_clean = str_remove_all(pep_clean, pattern = "\\.\\[.*\\]$")
pep_clean = sub("\\s+$", "", gsub('(.{5})', '\\1 ', pep_clean))
mod_sua$Clean_peptide = pep_clean

# order by histone, prot, and peptide position
start_peptide_pos = as.numeric(str_remove(str_split(mod_sua$Peptide_position, pattern = "-", simplify = T)[,1], "\\["))
mod_sua = mod_sua[order(mod_sua$Histone_Type, mod_sua$Protein_ID,start_peptide_pos ,mod_sua$Position_in_peptide, mod_sua$Modification),]

write.table(mod_sua, file = sprintf("%s.csv", out_fn), quote = F, sep = "\t", row.names = F)

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
