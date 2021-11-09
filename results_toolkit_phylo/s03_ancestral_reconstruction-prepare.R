# load libs
library(stringr)
library(ape)

# input
ort_fn = "orthogroups_euk.csv"
out_fn = "orthogroups_euk.ancestral"
phy_fn = "data/species_tree.newick"
count_path = "/home/xavi/Programes/Count/Count.jar"

# read phylogeny data
phy = read.tree(phy_fn)
sps_list = phy$tip.label


# #### PFAM-based training dataset ####
# input
dat_fn = "../data/pfam_annotations.csv.gz"

# read orthogroup data
dat = read.table(dat_fn, header = F, col.names = c("gene","orthogroup"))
dat$species = stringr::str_split(dat$gene, pattern = "_", simplify = T)[,1]
dat$species = factor(dat$species, levels = sps_list)

# output counts (all)
ort_counts = table(orthogroup = dat$orthogroup, species = dat$species)
write.table(ort_counts, file = "data/pfam_domain_counts.csv", sep="\t", quote = F,
            row.names = T, col.names = T)
ort_counts = read.table("data/pfam_domain_counts.csv")

# keep only domains present in at least 5% of the dataset
filt_rows = apply(ort_counts, 1, function(x) sum(x>0)) > ncol(ort_counts) * 0.05
ort_counts_filt = ort_counts[filt_rows,]

# get 1000 random rows
ort_counts_train = ort_counts_filt[sample(nrow(ort_counts_filt), 1000, replace = F),]
ort_counts_train[ort_counts_train>1] = 1

# output training set
ott_fn = "data/pfam_domain_counts.train.csv"
write.table(ort_counts_train, file = ott_fn, sep="\t", quote = F,
            row.names = T, col.names = T)
system(command = sprintf("sed -i 1's/^/\t/' %s", ott_fn))



#### POSSVM output ####

# read orthogroup data
ort = read.table(ort_fn, header = T)
ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
ort$species = factor(ort$species, levels = sps_list)

ort_counts = table(orthogroup = ort$orthogroup, species = ort$species)

# output counts (all)
write.table(ort_counts, file = sprintf("%s.counts.csv",out_fn), sep="\t", quote = F, 
            row.names = T, col.names = T)
system(command = sprintf("sed -i 1's/^/\t/' %s.counts.csv", out_fn))







stop("ara")

#### COUNT WAGNER ####

# run count
system(
  sprintf(
    "java -cp %s  ca.umontreal.iro.evolution.genecontent.AsymmetricWagner -gain 1 %s %s.counts.csv > %s.wagner_out.csv",
    count_path,
    phy_fn,
    out_fn,
    out_fn)
)

# get presence table
system(
  sprintf(
    "grep '^# FAMILY' %s.wagner_out.csv | cut -f 2-%i > %s.wagner_pres.csv",
    out_fn,
    length(phy$node.label)+length(phy$tip.label)+2,
    out_fn
  )
)

# get summary table
system(
  sprintf(
    "grep '^# FAMILY' %s.wagner_out.csv | cut -f 2,%i- > %s.wagner_summ.csv",
    out_fn,
    length(phy$node.label)+length(phy$tip.label)+3,
    out_fn
  )
)

