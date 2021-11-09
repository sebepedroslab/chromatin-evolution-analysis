# libraries
library(stringr)
library(seqinr)
library(msa)
library(Biostrings)

# load fusions
fas = Biostrings::readAAStringSet("fusions_all.fasta")
cla = read.table("histone_dimers_evidence.csv", sep = "\t", header = TRUE)

# load orthologs of zebrafish of a certain length
fos = Biostrings::readAAStringSet("fusions_orthologs-Zebrafish_si_ch211_113a14.fa")

# clean up a little bit: remove non-conserved genes and short proteins
black_list = c("ENSIPUP00000011897", "ENSGMOP00000010659", "ENSCARP00000018642")
fos = fos [ !names(fos) %in% black_list ]
fos = fos [ Biostrings::width(fos) > 200 ]
fos = unique(fos)



# concatentate
fas = c(fas,fos)


# define list of alignments to produce
list_ali = list(
  "met_H4-H2B" = cla [cla$fusion == "H4 H2B" , "gene"],
  "met_H2A-H4" = cla [cla$fusion == "H2A H4" , "gene"],
  "met_H3-H3" = cla [cla$fusion == "H3 H3" & cla$species != "Klenit" , "gene"],
  "ver_H2B-H2A" = c(cla [cla$species == "Drer" , "gene"], names(fos)),
  "amp_H2B-H2A" = cla [cla$species %in% c("Ddis","Ttra") , "gene"]
  
)

# align and draw
for (i in 1:length(list_ali)) {
  
  message("plot ", names(list_ali)[i])
  fai = fas [ list_ali[[i]] ]
  ali = msa::msa(fai, method = "ClustalW", maxiters = 1000, order = "input")
  
  # plot alignment
  msa::msaPrettyPrint(
    ali, askForOverwrite=FALSE, 
    shadingMode = "similar", 
    showConsensus = "none",
    showLogo = "none",
    shadingColors = "blues", 
    logoColors = "chemical",
    paperWidth = 32, paperHeight = 10, 
    output = "pdf",
    file = sprintf("aligned_fusions_%s.pdf", names(list_ali)[i]))
  
}
