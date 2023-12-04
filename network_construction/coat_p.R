# coat
source("/hwfssz5/ST_HEALTH/P18Z10200N0127/PROJECT/P18Z10200N0127_GRJ/juyanmei/bin/COAT/COAT/coat_v2.R")

args <- commandArgs(T)

prof <- read.csv(args[1], row.names = 1)
max_iter <- as.integer(args[2])
mini <- as.numeric(args[3])
prefix <- args[4]

# replace 0 and transpose prof
prof[prof == 0] <- mini
prof <- t(prof)

# coat
print("---------------Start coat------------------")
sim <- coat(prof, nFoler = 5, soft = 0.2)  #Is the adjacency matrix
colnames(sim) <- rownames(sim) <- colnames(prof)
write.csv(sim, paste(prefix, "/coat_Adj.csv", sep = ""))


#Bootstrap
print("Step 1: start bootstrap")
rows <- rownames(prof)
boot_arr <- array(dim = c(dim(sim)[1], dim(sim)[1], max_iter))
print(paste("iteration：", dim(boot_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  t_x <- sample(rows, replace = TRUE)
  sim <- coat(prof[t_x, ], nFoler = 10)
  boot_arr[, , i] <- sim
}
print("Bootstrap done")

#perm and renorm
print("Step 2: start permutation")
perm_arr <- array(dim = c(dim(sim)[1], dim(sim)[1], max_iter))
print(paste("iteration：", dim(perm_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  x <- apply(prof, 2, FUN = sample) #permutation
  x <- x / rowSums(x) #renormalization
  sim <- coat(x, nFoler = 10)
  perm_arr[, , i] <- sim
}
print("Permutation done")

#Converting NA to 0
boot_arr[is.na(boot_arr)] <- 0
perm_arr[is.na(perm_arr)] <- 0

p_val <- array(dim = c(dim(sim)[1], dim(sim)[1]))
for (i in 1:dim(sim)[1]) {
  for (j in 1:dim(sim)[1]) {
    t <- wilcox.test(perm_arr[i, j, ], boot_arr[i, j, ], alternative = "two.sided", paired = FALSE, exact = TRUE)
    p_val[i, j] <- t$p.value
  }
}

p <- p.adjust(p_val, method = "fdr")
dim(p) <- dim(p_val)
write.csv(p, paste(prefix, "/coat_Pval.csv", sep = ""))


print("-------------------Coat done-------------------------")