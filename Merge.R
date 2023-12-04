args <- commandArgs(T)

#library(gMCP)

scoring <- function(arr) {
  m <- max(abs(arr), na.rm = TRUE)
  arr <- (arr / m) * 100
  return(arr)
}

dir <- args[1]
cutoff <- args[2]
group <- args[3]

# merge score value
count <- 1
methods <- c("huge", "coat", "MI", "Bray-Curtis")
ori <- read.csv(paste(dir, "/", methods[1], "_Adj.csv", sep = ""), row.names = 1)
num <- nrow(ori)
m_score <- 0
p_m <- array(0, c(num, num, length(methods)))
w <- c(1/3, 1/3, 1/6, 1/6) #weights
for (i in methods) {
  data <- read.csv(paste(dir, "/", i, "_Adj.new.csv", sep = ""), row.names = 1)
  p <- read.csv(paste(dir, "/", i, "_Pval.new.csv", sep = ""), row.names = 1)
  p_m[, , count] <- as.matrix(p)
  m_score <- m_score + w[count] * abs(scoring(data))
  count <- count + 1
}

# Calculating the sign of the score
m_score_sign <- 0
methods <- c("huge", "coat")
count <- 1
w <- c(1/2, 1/2)
for (i in methods) {
  data <- read.csv(paste(dir, "/", i, "_Adj.new.csv", sep = ""), row.names = 1)
  m_score_sign <- m_score_sign + w[count] * scoring(data)
  count <- count + 1
}

m_score <- m_score * sign(m_score_sign)
rownames(m_score) <- colnames(m_score) <- rownames(ori)

#------------Merging p-values------------------
#library(gMCP)
# simes.test was copied from gMCP package which can not be installed.
simes.test <- function (pvalues, weights, alpha = 0.05, adjPValues = TRUE, verbose = FALSE, ...) {
    mJ <- Inf
    for (j in 1:length(pvalues)) {
        Jj <- pvalues <= pvalues[j]
        if (adjPValues) {
            mJt <- pvalues[j]/sum(weights[Jj])
            if (is.na(mJt)) {
                mJt <- 0
            }
            if (mJt < mJ) {
                mJ <- mJt
            }
        }
    }
    if (adjPValues) {
        return(mJ)
    }
    else {
        return(mJ <= alpha)
    }
}

sime <- function(p) {
  x <- simes.test(p, weights = c(1/3, 1/3, 1/6, 1/6))
  return(x)
}

p_f <- apply(p_m, c(1, 2), sime)
#p_f[p_f == 0] <- NA
#p-values with NA become 0 because if NA is present it is assigned 0
rownames(p_f) <- colnames(p_f) <- rownames(ori)

#write.csv(p_f, paste(group, "_pvalues.csv", sep = ""))

score_pval_filt <- m_score
score_pval_filt[p_f > cutoff] <- 0
write.csv(score_pval_filt, paste(group, "_merge_scores.csv", sep = ""))
