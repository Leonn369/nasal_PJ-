args <- commandArgs(T)

library(huge)
library(mboost)
library(boot)
library(reshape2)

prof <- read.csv(args[1], row.names = 1)
max_iter <- as.integer(args[2])
prefix <- args[3]

# Huge
print("--------------Start Huge----------------")
prof <- as.data.frame(t(prof))
prof_trans <- huge.npn(prof)
res <- huge(prof_trans, nlambda = 100)
select <- huge.select(res)
res_opt <- huge(prof_trans, lambda = select$opt.lambda)
res_opt_beta <- as.matrix(res_opt$beta[[1]])
colnames(res_opt_beta) <- rownames(res_opt_beta) <- colnames(prof)

#Bootstrap distribution
#------------Using boot-----------------
print("Step 1: start bootstrap")
set.seed(0)
bootStat <- function(data, indices) {
  data <- data[indices, ]
  data_trans <- huge.npn(data)
  res <- huge(data_trans, nlambda = 100)
  select <- huge.select(res)
  a <- huge(data_trans, lambda = select$opt.lambda)
  b <- matrix(a$beta[[1]])
  return(b)
}
model_boot <- boot(prof, bootStat, max_iter)
print(paste("iteration：", max_iter, sep = ""))
print("Bootstrap done")

#-----------Permutation with renormalization-------------
print("Step 2: start permutation")
#copy of the data
prof_p <- as.data.frame(prof)
out_comb <- matrix(0, dim(model_boot$t)[2], 1)
out_comb <- NULL
#permutation
counter <- 0
set.seed(0)
while (counter < max_iter) {
  prof_p <- sample(prof_p)
  prof_p <- (prof_p / rowSums(prof_p)) #renormalization
  out <- bootStat(prof_p, indices = 1:dim(prof)[1]) # nolint
  out_comb <- cbind(out_comb, melt(as.matrix(out))$value)
  counter <- counter + 1
}
print(paste("iteration：", max_iter, sep = ""))
print("Permutation done")

#Comparing two distributions
p_test <- c()
for (i in 1:dim(out_comb)[1]) {
  p <- wilcox.test(model_boot$t[, i], out_comb[i, ], alternative = "two.sided", paired = FALSE)$p.value # nolint
  p_test <- c(p_test, p)
}
#correction of multiple comparision
p_test <- p.adjust(p_test, method = "fdr")
p_m <- matrix(p_test, ncol(prof), ncol(prof))
colnames(p_m) <- rownames(p_m) <- colnames(data)
write.csv(res_opt_beta, paste(prefix, "/huge_Adj.csv", sep = ""))
write.csv(p_m, paste(prefix, "/huge_Pval.csv", sep = ""))


print("------------------- Huge done-------------------------")
