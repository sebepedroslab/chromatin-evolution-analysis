# functions and libraries
library(stringr)
library(seqinr)
library(readxl)


# load coverage
list_xlsxs = list.files(path = "../../results_PTMs/data-archaea/", pattern = "*me-ac.xlsx", full.names = TRUE)
# list_xlsxs = list.files(path = "../../results_PTMs/data-archaea/", pattern = "*.xlsx", full.names = TRUE)
for (i in 1:length(list_xlsxs)) {
    
    fn = list_xlsxs[i]
    dat = read_xlsx(fn, sheet = 1)
    daf = dat [ grep("Histone", dat$Accession) ,]
    
    if (i == 1) { 
        dao = daf   
    } else {
        dao = dao [ , colnames(daf) ]
        dao = rbind(dao, daf)
    }
    
}

# load sequences
list_fasta = list.files(path = "proteomes/", pattern = "*.fasta", full.names = TRUE)
for (i in 1:length(list_fasta)) {
    
    fn = list_fasta[i]
    if (i == 1) {
        fas = seqinr::read.fasta(fn)
    } else {
        fas = c(fas, seqinr::read.fasta(fn))
    }
    
}

fao = fas [ which(names(fas) %in% dao$Accession) ]
fao = fao [ dao$Accession ]
fao_i = lapply(fao, function(i) toupper(paste(i, collapse = "")))
fao_f = lapply(fao_i, function(i) str_count(i, "K") / str_length(i) )

pdf("coverage_pi_ph.pdf", width = 10, height = 2.5)
par(mar=c(5,12,5,2))
layout(mat = matrix(1:4, nrow = 1))

# length
b=barplot(rev(dao$`# AAs`), horiz = TRUE, las = 1,    xlim = c(0,100), main = "Protein length", xlab = "# aa",
        names.arg = rev(dao$Accession), cex.names = 0.5)
text(x=0, y=b, sprintf(rev(dao$`# AAs`), fmt = "%i"), col = "red", cex = 0.7, pos = 4)

# coverage
b=barplot(rev(dao$`Coverage [%]`), horiz = TRUE, las = 1, xlim = c(0,100), main = "Coverage %", xlab = "%")
text(x=0, y=b, sprintf(rev(dao$`Coverage [%]`), fmt = "%i"), col = "red", cex = 0.7, pos = 4)

# fraction Ks
b=barplot(rev(unlist(fao_f)), horiz = TRUE, las = 1, xlim = c(0,0.2),  main = "% K residues", xlab = "%", names.arg = NA)
text(x=0, y=b, sprintf(rev(fao_f), fmt = "%.3f"), col = "red", cex = 0.7, pos = 4)

# isoelectric point
b=barplot(rev(dao$`calc. pI`), horiz = TRUE, las = 1, xlim = c(0,12),  main = "isoelectric point", xlab = "pI")
text(x=0, y=b, sprintf(rev(dao$`calc. pI`), fmt = "%.2f"), col = "red", cex = 0.7, pos = 4)

par(mar=c(5,2,5,2))
plot(x=rep(NA, length(b)), y=b, xlim = c(0,1), frame.plot=FALSE, axes = FALSE)
text(x=0, y=b, paste(rev(dao$Accession), rev(dao$Modifications)), cex = 0.4, pos = 4)
dev.off()


# save
write.table(dao, "coverage_pi_ph.csv", sep = "\t", row.names = FALSE, quote = FALSE)
