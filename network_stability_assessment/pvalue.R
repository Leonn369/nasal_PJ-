args <- commandArgs(T)

library(igraph)
library(pulsar)
library(influential)
library(reshape2)
library(stringi)
library(stringr)
source("swan.R")
library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

cut_cor <- function(data, cutoff) {
  data_cut <- data
  data_cut[abs(data_cut) < cutoff] <- 0
  return(data_cut)
}

ivi_random_cluster <- function(data,  k = 100, d = d,  norm = T, total = T, positive = F) {
  colnames(data) <- rownames(data) <- str_replace(colnames(data), 'd__Fungi.+.g__', 'd__Fungi|g__')
  mm <- as.matrix(data)
  if (total) {
    mm <- abs(mm)
  } else {
    if (positive) {
      mm[mm < 0] <- 0
    } else {
      mm[mm > 0] <- 0
      mm <- abs(mm)
    }
  }
  mm2 <- mm
  mm2[mm!=0] <- 1
  Mygraph <- graph_from_adjacency_matrix(mm2, mode = 'undirected', weighted = T)
  print(Mygraph)
  f1 <- swan_combinatory_weight(mm, d = d, k = k, norm = T)
  return(f1)
}
cluster <- as.character(args[1])
cutoff <- as.numeric(args[2])
d <- as.numeric(args[3])
total <- args[4]
positive <- args[5]
prefix <- args[6]
res <- args[7]

# original data
if (cluster == "male") {
  cyto_cluster1 <- read.csv(paste(prefix, 'male_merge_scores.csv', sep = ''), header = T, row.names = 1, check.names = F)
  cyto_cluster1_cut <- cut_cor(cyto_cluster1, cutoff)/100
  cluster1_perm <- ivi_random_cluster(cyto_cluster1_cut, d = d, k = 1000, norm = F, total = total, positive = positive)
  write.csv(cluster1_perm$fin, paste(res, "/fin_perm_cluster1_cutoff", cutoff, "_d", d, ".csv", sep = ""))
  write.csv(cluster1_perm$fin_p, paste(res, "/fin_pval_perm_cluster1_cutoff", cutoff, "_d", d, ".csv", sep = ""))
  write.csv(cluster1_perm$rand, paste(res, "/random_perm_cluster1_cutoff", cutoff, "_d", d, ".csv", sep = ""))
  write.csv(cluster1_perm$fin_statistic, paste(res, "/fin_statistic_perm_cluster1_cutoff", cutoff, "_d", d, ".csv", sep = ""))
} else if (cluster == "female") {
  cyto_cluster2 <- read.csv(paste(prefix, 'female_merge_scores.csv', sep = ''), header = T, row.names = 1, check.names = F) 
  cyto_cluster2_cut <- cut_cor(cyto_cluster2, cutoff)/100
  cluster2_perm <- ivi_random_cluster(cyto_cluster2_cut, d = d, k = 1000, norm = F, total = total, positive = positive)
  write.csv(cluster2_perm$fin, paste(res, "/fin_perm_cluster2_cutoff", cutoff, "_d", d, ".csv", sep = ""))  
  write.csv(cluster2_perm$fin_p, paste(res, "/fin_pval_perm_cluster2_cutoff", cutoff, "_d", d, ".csv", sep = ""))
  write.csv(cluster2_perm$rand, paste(res, "/random_perm_cluster2_cutoff", cutoff, "_d", d, ".csv", sep = ""))
  write.csv(cluster2_perm$fin_statistic, paste(res, "/fin_statistic_perm_cluster2_cutoff", cutoff, "_d", d, ".csv", sep = ""))
}

stopCluster(cl)
