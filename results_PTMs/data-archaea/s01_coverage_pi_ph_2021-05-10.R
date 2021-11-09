# functions and libraries
library(stringr)
library(seqinr)
library(readxl)


# load coverage
# list_xlsxs = list.files(path = ".", pattern = "*me-ac.xlsx", full.names = TRUE)
list_xlsxs = list.files(path = ".", pattern = "*.xlsx", full.names = TRUE)
for (i in 1:length(list_xlsxs)) {
    
    fn = list_xlsxs[i]
    dat = read_xlsx(fn, sheet = 2)
    daf = dat [ grep("Histone", dat$`Master Protein Accessions`) ,]
    daf = daf [ !is.na(daf$`Modifications in Master Proteins`) ,]
    daf = daf [ grep("Acetyl|Methyl",daf$`Modifications in Master Proteins`) ,]
    
    if (i == 1) { 
        dao = daf   
    } else {
        dao = dao [ , colnames(daf) ]
        dao = rbind(dao, daf)
    }
    
}


write.table(dao, "../modifications_per_sps_archaea.csv", row.names = FALSE, quote = FALSE, sep = "\t")