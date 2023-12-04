args <- commandArgs(T)

statstic_cor <- function(data_r, data_p, cut = 0.01) {
  r <- data_r
  p <- data_p
  p[is.na(p)] <- 1
  r[is.na(r)] <- 0
  for (i in 1:nrow(r)) {
    for (j in 1:ncol(r)) {
      if (p[i, j] > cut) {
        r[i, j] <- 0
      }
    }
  }
  cor_ul <- data.frame(upper = r[upper.tri(r)],
                    lower = t(r)[upper.tri(r)])
  p_ul <- data.frame(upper = p[upper.tri(p)],
                    lower = t(p)[upper.tri(p)])
  res <- c(0, 0, 0, 0, 0)
  names(res) <- c("equal0", "equalNot0", "same_direction", "one_0", "differ_direction")
  for (i in 1:nrow(cor_ul)) {
    cor_ul[i, 1] <- round(cor_ul[i, 1], 8)
    cor_ul[i, 2] <- round(cor_ul[i, 2], 8)
    if (cor_ul[i, 1] == cor_ul[i, 2]) {
      if (cor_ul[i, 1] == 0) {
        res[1] <- res[1] + 1
      } else {
        res[2] <- res[2] + 1
      }
    } else {
      if (cor_ul[i, 1] * cor_ul[i, 2] > 0) {
        res[3] <- res[3] + 1
      } else if (cor_ul[i, 1] * cor_ul[i, 2] == 0) {
         res[4] <- res[4] + 1
      } else if (cor_ul[i, 1] * cor_ul[i, 2] < 0) {
         res[5] <- res[5] + 1
      }
    }
  }
  cor_ul_mod <- matrix(0, ncol(data_r), ncol(data_r))
  cor_ul_mod[upper.tri(cor_ul_mod)] <- apply(cor_ul, 1, mean)
  cor_ul_mod[lower.tri(cor_ul_mod)] <- t(cor_ul_mod)[lower.tri(cor_ul_mod)]
  
  p_ul_mod <- matrix(1, ncol(data_r), ncol(data_r))
  p_ul_mod[upper.tri(p_ul_mod)] <- apply(p_ul, 1, mean)
  p_ul_mod[lower.tri(p_ul_mod)] <- t(p_ul_mod)[lower.tri(p_ul_mod)]

  final <- list(est = cor_ul_mod,
                pvalue = p_ul_mod,
                ratio = res / (ncol(r) * (ncol(r) - 1)) * 100 * 2)
  return(final)
}

dir <- args[1]
pvalue <- as.numeric(args[2])
prefix <- args[3]
group <- args[4]

methods <- c("coat", "huge", "MI", "Bray-Curtis")
methods_n <- length(methods)
num <- nrow(read.csv(paste(dir, "/", methods[1], "_Adj.csv", sep = ""), row.names = 1))
p_m <- array(0, c(num, num, methods_n))
pp_m <- array(1, c(num, num, methods_n))
count <- 1
ratio <- NULL

for (i in methods) {
  data <- read.csv(paste(dir, "/", i, "_Adj.csv", sep = ""), row.names = 1)
  p <- read.csv(paste(dir, "/", i, "_Pval.csv", sep = ""), row.names = 1)
  print(paste(i, "Start", sep = " "))
  k <- statstic_cor(data, p, cut = pvalue)
  p_m[, , count] <- k$est
  pp_m[, , count] <- k$pvalue
  ratio <- rbind(ratio, k$ratio)
  rownames(p_m[, , count]) <- colnames(p_m[, , count]) <- rownames(data)
  rownames(pp_m[, , count]) <- colnames(pp_m[, , count]) <- rownames(data)
  write.csv(p_m[, , count], paste(dir, "/", i, "_Adj.new.csv", sep = "")) # nolint
  write.csv(pp_m[, , count], paste(dir, "/", i, "_Pval.new.csv", sep = ""))  # nolint
  count <- count + 1
  print(paste(i, "Done", sep = " "))
}
rownames(ratio) <- methods
write.csv(ratio, paste(prefix, "/", group, "_ratio.csv", sep = ""))

pdf(paste(prefix, "/", group, "_integration_hist.pdf", sep = ""))
par(mfrow = c(2, 2))
for (i in 1:length(methods)) {
  p_m_max <- round(max(abs(p_m[, , i])), 3)
  p_m_min <- round(min(abs(p_m[, , i])), 3)
  p_m_median <- round(median(abs(p_m[, , i])), 3)
  hist(p_m[, , i], xlab = paste("max:", p_m_max, "; min:", p_m_min, "; median:", p_m_median, sep = ""), main = methods[i])
}
dev.off()
