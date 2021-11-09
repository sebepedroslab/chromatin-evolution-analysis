# libraries
library(stringr)
library(seqinr)
library(msa)
library(Biostrings)

# define list of alignments to produce
list_ali = list(
  "stamsii_modified.g.fasta"
)

# align and draw
for (i in 1:length(list_ali)) {
  
  fai = list_ali[[i]]
  ali = Biostrings::readAAMultipleAlignment(fai, format = "fasta")
  
  # plot alignment
  msa::msaPrettyPrint(
    ali, askForOverwrite=FALSE, 
    shadingMode = "functional", 
    shadingModeArg = "chemical",
    showConsensus = "none",
    showLogo = "none",
    shadingColors = "blues", 
    logoColors = "chemical",
    psFonts = TRUE,
    paperWidth = 32, paperHeight = 10, 
    output = "pdf",
    file = sprintf("%s.colour.pdf",basename(fai)))
  
}
