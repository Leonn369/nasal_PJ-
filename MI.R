args <- commandArgs(T)
library("minet")

prof <- read.csv(args[1], row.names = 1)
max_iter <- as.integer(args[2])
prefix <- args[3]


#Mutual Information
prof <- as.data.frame(t(prof))
print("--------------Start MI----------------")
m <- build.mim(prof, estimator = "mi.mm", disc = "equalwidth")
net <- aracne(m)
write.csv(net, paste(prefix, "/MI_Adj.csv", sep = ""))

#Bootstrap
print("Step 1: start bootstrap")
rows <- rownames(prof)
boot_arr <- array(dim = c(dim(net)[1], dim(net)[1], max_iter))
print(paste("iteration：", dim(boot_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  t_x <- sample(rows, replace = TRUE)
  m <- build.mim(prof[t_x, ], estimator = "mi.mm", disc = "equalwidth")
  net <- aracne(m)
  boot_arr[, , i] <- net
}
print("Bootstrap done")

#perm and renorm
print("Step 2: start permutation")
perm_arr <- array(dim = c(dim(net)[1], dim(net)[1], max_iter))
print(paste("iteration：", dim(perm_arr)[3], sep = ""))
set.seed(0)
for (i in 1:max_iter) {
  x <- apply(prof, 2, FUN = sample) #permutation
  x <- x / rowSums(x) #renormalization
  m <- build.mim(x, estimator = "mi.mm", disc = "equalwidth")
  net <- aracne(m)
  #print(sum(is.na(net)))
  perm_arr[, , i] <- net
}
print("Permutation done")

#Converting NA to 0
boot_arr[is.na(boot_arr)] <- 0
perm_arr[is.na(perm_arr)] <- 0

p_val = array(dim = c(dim(net)[1], dim(net)[1]))
for (i in 1:dim(net)[1]) {
  for (j in 1:dim(net)[1]) {
    t <- wilcox.test(perm_arr[i, j, ], boot_arr[i, j, ], alternative = "two.sided", paired = FALSE, exact = TRUE)
    p_val[i, j] <- t$p.value
  }
}

p <- p.adjust(p_val, method = "fdr")
dim(p) <- dim(p_val)
write.csv(p, paste(prefix, "/MI_Pval.csv", sep = ""))


print("------------------- MI done-------------------------")