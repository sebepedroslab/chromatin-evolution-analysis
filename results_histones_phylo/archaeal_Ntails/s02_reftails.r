# libraries
library(stringr)
library(seqinr)
library(data.table)



pdf("dist_euk_histones_tail_Kfreq.pdf", width = 4, height = 4)

# histone 3
chs = read.fasta("canonical_H3.fasta")
chs_i = chs
chs_i = lapply(chs_i, function(i) toupper(paste(i, collapse = "")))
chs_i = lapply(chs_i, function(i) gsub("-","",i))
chs_i = lapply(chs_i, function(i) gsub("PGTVAL.*", "PGT", i))
chs_i_len = lapply(chs_i, function(i) str_length(i))
chs_i_nks = lapply(chs_i, function(i) str_count(i, "K") )

# plot
vioplot::vioplot(unlist(chs_i_nks)/unlist(chs_i_len), col = "blue", ylim = c(0.1,0.3))
title(sub = sprintf("median f = %.3f | l = %i", median(unlist(chs_i_nks)/unlist(chs_i_len)), median(unlist(chs_i_len))), main = "H3")
abline(h = median(unlist(chs_i_nks)/unlist(chs_i_len)), lty = 2)

# histone 4
chs = read.fasta("canonical_H4.fasta")
chs_i = chs
chs_i = lapply(chs_i, function(i) toupper(paste(i, collapse = "")))
chs_i = lapply(chs_i, function(i) gsub("-","",i))
chs_i = lapply(chs_i, function(i) gsub("K[VIS][LN][RSK][DA].*", "", i))
chs_i_len = lapply(chs_i, function(i) str_length(i))
chs_i_nks = lapply(chs_i, function(i) str_count(i, "K") )

# plot
vioplot::vioplot(unlist(chs_i_nks)/unlist(chs_i_len), col = "blue", ylim = c(0.1,0.3))
title(sub = sprintf("median f = %.3f | l = %i", median(unlist(chs_i_nks)/unlist(chs_i_len)), median(unlist(chs_i_len))), main = "H4")
abline(h = median(unlist(chs_i_nks)/unlist(chs_i_len)), lty = 2)

dev.off()


