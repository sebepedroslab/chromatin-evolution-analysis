# load libs
library(Rphylip)
library(stringr)
library(ape)

# input
phy_fn = "data/species_tree.newick"
ort_fn = "orthogroups_euk.csv"

# read phylogeny data
phy = read.tree(phy_fn)
sps_list = phy$tip.label
phy$edge.length = 1


# read orthogroup data
ort = read.table(ort_fn, sep="\t", header = T)
ort$species = stringr::str_split(ort$gene, pattern = "_", simplify = T)[,1]
ort$species = factor(ort$species, levels = sps_list)
ort$orthogroup_id = stringr::str_split(ort$orthogroup, pattern = ":", simplify = T)[,1]
ort_list = unique(ort$orthogroup_id)
ort_counts = as.data.frame(table(species = ort$species, orthogroup_id = ort$orthogroup_id))

# get one vector
ort_ix = ort_list == "Hist_deacetyl.HG2.6:like:HDAC1/HDAC2/HDAC3"
ori = ort_counts[ort_counts$orthogroup_id == "Hist_deacetyl.HG2.6", "Freq" ]
ori = as.numeric(ori>0)
names(ori) = sps_list
# matrix
orm = as.matrix(ori, ncol=1)
rownames(orm) = sps_list


orm_tre = Rphylip::Rpars(X=orm, path = "/home/xavi/Programes/miniconda3/bin/", cleanup=T)
orm_par = Rphylip::Rpars(X=orm, path = "/home/xavi/Programes/miniconda3/bin/", cleanup=T, tree=orm_tre)

stop("THIS")

# load phytools
library(phytools)

# input
phy_fn = "data/species_tree.newick"
ort_fn = "orthogroups_euk.csv"

# read phylogeny data
phy = read.tree(phy_fn)
sps_list = phy$tip.label




# get transition matrix
Q = matrix(c(0, 1,  10, 0), 2, 2)
rownames(Q) = c("0","1")
colnames(Q) = c("0","1")
Q

# phyfitER = rerootingMethod_mod(phy, ori, model = "SYM")
fitER = rerootingMethod(phy, as.matrix(ori), model = "ER")


fitMk = fitMk(phy, ori,model="ER")
plot.fitMk(fitMk)


## simulate a tree & some data
tree <- pbtree(n = 26, scale = 1)
tree$tip.label <- LETTERS[26:1]
Q <- matrix(c(-2, 1, 1, 1, -2, 1, 1, 1, -2), 3, 3)
colnames(Q) <- rownames(Q) <- letters[1:3]
x <- sim.history(tree, Q, anc = "a")$states
print(x)

fitER <- rerootingMethod(tree, x, model = "ER")
print(fitER)






# function

rerootingMethod_mod = function (tree, x, model = c("ER", "SYM"))  {
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (hasArg(tips)) 
    tips <- list(...)$tips
  else tips <- NULL
  if (!is.matrix(model)) 
    model <- model[1]
  n <- Ntip(tree)
  if (!is.matrix(x)) {
    yy <- to.matrix(x, sort(unique(x)))
    if (is.null(tips)) 
      tips <- FALSE
  }
  else {
    if (is.null(tips)) 
      tips <- TRUE
    yy <- x
  }
  yy <- yy[tree$tip.label]
  YY <- fitMk(tree, yy, model = model, output.liks = TRUE)
  Q <- matrix(c(0, YY$rates)[YY$index.matrix + 1], length(YY$states), 
              length(YY$states), dimnames = list(YY$states, YY$states))
  diag(Q) <- -colSums(Q, na.rm = TRUE)
  nn <- if (tips) 
    c(1:n, if (tree$Nnode > 1) 2:tree$Nnode + n)
  else {
    if (tree$Nnode > 1) 
      2:tree$Nnode + n
    else vector()
  }
  ff <- function(nn) {
    tt <- if (nn > Ntip(tree)) 
      ape::root.phylo(tree, node = nn)
    else
      reroot(tree, nn, tree$edge.length[which(tree$edge[, 2] == nn)])
    fitMk(tt, yy, model = model, fixedQ = Q, output.liks = TRUE)$lik.anc[1, 
    ]
  }
  XX <- t(sapply(nn, ff))
  if (tips) 
    XX <- rbind(XX[1:n, ], YY$lik.anc[1, ], if (tree$Nnode > 
                                                1) 
      XX[(n + 1):nrow(XX), ])
  else XX <- rbind(yy, YY$lik.anc[1, ], if (tree$Nnode > 1) 
    XX)
  rownames(XX) <- 1:(tree$Nnode + n)
  if (tips) 
    rownames(XX)[1:n] <- tree$tip.label
  XX <- if (tips) 
    XX
  else XX[1:tree$Nnode + n, ]
  obj <- list(loglik = YY$logLik, Q = Q, marginal.anc = XX, 
              tree = tree, x = yy)
  class(obj) <- "rerootingMethod"
  obj
}
