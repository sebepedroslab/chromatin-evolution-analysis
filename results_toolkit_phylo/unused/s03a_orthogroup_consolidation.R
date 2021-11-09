library(stringr)

# load orthogroups
ogs_fn = "orthogroups_euk.csv"
ogs = read.table(ogs_fn, header = T, sep="\t", stringsAsFactors = F)


# get affinity dict: table relating each OG with a "like" particle to its reference OG
ogs_affinity = data.frame(
  orthogroup = ogs[!duplicated(ogs$orthogroup) & grepl(pattern=":likeclu:", x=ogs$orthogroup),]$orthogroup,
  stringsAsFactors = F
)

# get original clu name
ogs_affinity$is_like_clu = stringr::str_split(ogs_affinity$orthogroup, pattern =  ":", simplify = T)[,5]
ogs_affinity$is_like_gen = stringr::str_split(ogs_affinity$orthogroup, pattern =  ":", simplify = T)[,3]
ogs_affinity$is_like_ogi = paste(
  gsub(pattern="\\.\\d+",replacement = "", x=stringr::str_split(ogs_affinity$orthogroup, pattern =  ":", simplify = T)[,1]),
  ".",ogs_affinity$is_like_clu,
  ":",ogs_affinity$is_like_gen,
  sep=""
)
# discard clusters that resemble more than one reference cluster
ogs_affinity = ogs_affinity[!grepl(pattern="/", x=ogs_affinity$is_like_clu),]


# replacements
ogs_new = ogs
ogs_new$orthogroup = apply(ogs, 1, FUN = function(x) {
  r = ogs_affinity[ogs_affinity$orthogroup == x[2], "is_like_ogi"]
  if (length(r)>0){ 
    x[2] = r
  } else {
    x[2] = x[2]
  }
})


# store new orthogroups
write.table(ogs_new, file = "orthogroups_euk.consolidated.csv", sep="\t", quote = F, row.names = F)

