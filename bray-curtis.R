# Bray-curtis
library(vegan)

args <- commandArgs(T)

prof <- read.csv(args[1], row.names = 1)
max_iter <- as.integer(args[2])
prefix <- args[3]

#Bray-Curtis
print("--------------Start Bray-Curtis----------------")
dis <- vegdist(prof, method = "bray", diag = TRUE, upper = TRUE)
sim <- (as.matrix(dis) - 1) * -1  #Is the adjacency matrix
write.csv(sim, paste(prefix, "/Bray-Curtis_Adj.csv", sep = ""))


#Bootstrap
print("Step 1: start bootstrap")
cols <- colnames(prof)
boot_arr <- array(dim = c(dim(sim)[1], dim(sim)[1], max_iter))
print(paste("iteration：", dim(boot_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  t_x <- sample(cols, replace = TRUE)
  dis <- vegdist(prof[,t_x], method = "bray", diag = TRUE, upper = TRUE)
  sim <- (as.matrix(dis) - 1) * -1
  boot_arr[, , i] <- sim
}
print("Bootstrap done")

#perm and renorm
print("Step 2: start permutation")
perm_arr <- array(dim = c(dim(sim)[1], dim(sim)[1], max_iter))
print(paste("iteration：", dim(perm_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  x <- t(apply(prof, 1, FUN = sample)) #permutation
  x <- t(t(x) / colSums(x)) #renormalization
  dis <- vegdist(x, method = "bray", diag = TRUE, upper = TRUE)
  sim <- (as.matrix(dis) - 1) * -1
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
write.csv(p, paste(prefix, "/Bray-Curtis_Pval.csv", sep = ""))

print("-------------------Bray-Curtis done-------------------------")