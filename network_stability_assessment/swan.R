swan_combinatory_unweight <- function(g_dat, k, d, norm = T) {
  g <- graph_from_adjacency_matrix(g_dat, mode = 'undirected', weighted = NULL)
  n <- length(V(g))
  noeud <- V(g)
  arc <- E(g)
  fin <- matrix(ncol = 6, nrow = n+1, 0)
  fin_p <- matrix(ncol = 6, nrow = n+1, 1)
  for(i in 1:n) {
    fin[i, 1] <- (i-1)/n
    fin_p[i, 1] <- (i-1)/n
  }
  fin[n+1, 1] <- 1
  fin_p[n+1, 1] <- 1
  
  # random
  test_function_normT <- function(p) {
    library(igraph)
    library(pulsar)
    res <- matrix(0, n, 1)
    al <- sample(1:n, n)
    g2 <- g
    for (i in 1:n) {
      g2_dat <- as_adjacency_matrix(g2)
      res[i, 1] <- natural.connectivity(g2_dat, norm = T)
      g2 <- delete_vertices(g2, al[i])
      al[al > al[i]] <- al[al > al[i]] - 1
    }
    return(res)
  }
  test_function_normF <- function(p) {
    library(igraph)
    library(pulsar)
    res <- matrix(0, n, 1)
    al <- sample(1:n, n)
    g2 <- g
    for (i in 1:n) {
      g2_dat <- as_adjacency_matrix(g2)
      res[i, 1] <- natural.connectivity(g2_dat, norm = F)
      g2 <- delete_vertices(g2, al[i])
      al[al > al[i]] <- al[al > al[i]] - 1
    }
    return(res)
  }
  if(norm) {
    randomMatrix <- matrix(unlist(foreach(x = c(1:k)) %dopar% test_function_normT(x)), n, k)
    fin[, 6] <- c(apply(randomMatrix, 1, mean), 0)
  } else {
    randomMatrix <- matrix(unlist(foreach(x = c(1:k)) %dopar% test_function_normF(x)), n, k)
    fin[, 6] <- c(apply(randomMatrix, 1, mean), 0)
  }
  # betweeness
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  bet <- betweenness(g)
  mat[, 2] <- bet
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2)
    fin[i, 2] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > 
                                                  matri[v, 1], 1] - 1
  }
  # degree
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  deg <- degree(g)
  mat[, 2] <- deg
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2)
    fin[i, 3] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > 
                                                  matri[v, 1], 1] - 1
  }
  # closeness
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  clo <- closeness(g)
  mat[, 2] <- clo
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2)
    fin[i, 4] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > 
                                                  matri[v, 1], 1] - 1
  }
  # ivi
  g2 <- g
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  ivi_data <- ivi(g2, vertices = noeud, directed = FALSE, d = d,  scaled = TRUE)
  mat[, 2] <- ivi_data
  matri <- mat[order(mat[, 2]), ]
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2)
    fin[i, 5] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > 
                                                  matri[v, 1], 1] - 1
  }
  for(i in 1:n) {
    fin_p[i, 2] <- length(which(randomMatrix[i, ] > fin[i, 2]))/k
    fin_p[i, 3] <- length(which(randomMatrix[i, ] > fin[i, 3]))/k
    fin_p[i, 4] <- length(which(randomMatrix[i, ] > fin[i, 4]))/k
    fin_p[i, 5] <- length(which(randomMatrix[i, ] > fin[i, 5]))/k
    fin_p[i, 6] <- length(which(randomMatrix[i, ] > fin[i, 6]))/k
  }
  fin <- as.data.frame(fin)
  colnames(fin) <- c("node_ratio", "betweenness", "degree", "closeness", "ivi", "random")
  colnames(fin_p) <- colnames(fin)
  lst <- list(NULL)
  lst$fin <- fin
  lst$fin_p <- fin_p
  lst$rand <- randomMatrix
  return(lst)
}

swan_combinatory_weight <- function(g_dat, k, d, norm = T) {
  g <- graph_from_adjacency_matrix(g_dat, mode = 'undirected', weighted = T)
  n <- length(V(g))
  noeud <- V(g)
  arc <- E(g)
  fin <- matrix(ncol = 6, nrow = n+1, 0)
  fin_p <- matrix(ncol = 6, nrow = n+1, 1)
  fin_statistic <- matrix(ncol = 4, nrow = n, 1)
  for(i in 1:n) {
    fin[i, 1] <- (i-1)/n
    fin_p[i, 1] <- (i-1)/n
  }
  fin[n+1, 1] <- 1
  fin_p[n+1, 1] <- 1
  
  # random
  test_function_normT <- function(p) {
    library(igraph)
    library(pulsar)
    res <- matrix(0, n, 1)
    al <- sample(1:n, n)
    g2 <- g
    for (i in 1:n) {
      g2_dat <- as_adjacency_matrix(g2, attr = "weight")
      res[i, 1] <- natural.connectivity(g2_dat, norm = T)
      g2 <- delete_vertices(g2, al[i])
      al[al > al[i]] <- al[al > al[i]] - 1
    }
    return(res)
  }
  test_function_normF <- function(p) {
    library(igraph)
    library(pulsar)
    res <- matrix(0, n, 1)
    al <- sample(1:n, n)
    g2 <- g
    for (i in 1:n) {
      g2_dat <- as_adjacency_matrix(g2, attr = "weight")
      res[i, 1] <- natural.connectivity(g2_dat, norm = F)
      g2 <- delete_vertices(g2, al[i])
      al[al > al[i]] <- al[al > al[i]] - 1
    }
    return(res)
  }
  if(norm) {
    randomMatrix <- matrix(unlist(foreach(x = c(1:k)) %dopar% test_function_normT(x)), n, k)
    fin[, 6] <- c(apply(randomMatrix, 1, mean), 0)
  } else {
    randomMatrix <- matrix(unlist(foreach(x = c(1:k)) %dopar% test_function_normF(x)), n, k)
    fin[, 6] <- c(apply(randomMatrix, 1, mean), 0)
  }
  
  # betweeness
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  bet <- betweenness(g)
  mat[, 2] <- bet
  fin_statistic[, 1] <- bet
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2, attr = "weight")
    fin[i, 2] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > matri[v, 1], 1] - 1
  }
  # degree
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  deg <- degree(g)
  mat[, 2] <- deg
  fin_statistic[, 2] <- deg
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2, attr = "weight")
    fin[i, 3] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > matri[v, 1], 1] - 1
  }
  # closeness
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  clo <- closeness(g)
  mat[, 2] <- clo
  fin_statistic[, 3] <- clo
  matri <- mat[order(mat[, 2]), ]
  g2 <- g
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2, attr = "weight")
    fin[i, 4] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > matri[v, 1], 1] - 1
  }
  # ivi
  g2 <- g
  mat <- matrix(ncol = 2, nrow = n, 0)
  mat[, 1] <- 1:n
  ivi_data <- ivi(g2, vertices = noeud, directed = FALSE, d = d,  scaled = TRUE)
  mat[, 2] <- ivi_data
  fin_statistic[, 4] <- ivi_data
  matri <- mat[order(mat[, 2]), ]
  for (i in 1:n) {
    v = n + 1 - i
    g2_dat <- as_adjacency_matrix(g2, attr = "weight")
    fin[i, 5] <- natural.connectivity(g2_dat, norm = norm)
    g2 <- delete_vertices(g2, matri[v, 1])
    matri[matri[, 1] > matri[v, 1], 1] <- matri[matri[, 1] > matri[v, 1], 1] - 1
  }
  for(i in 1:n) {
    fin_p[i, 2] <- length(which(randomMatrix[i, ] > fin[i, 2]))/k
    fin_p[i, 3] <- length(which(randomMatrix[i, ] > fin[i, 3]))/k
    fin_p[i, 4] <- length(which(randomMatrix[i, ] > fin[i, 4]))/k
    fin_p[i, 5] <- length(which(randomMatrix[i, ] > fin[i, 5]))/k
    fin_p[i, 6] <- length(which(randomMatrix[i, ] > fin[i, 6]))/k
  }
  fin <- as.data.frame(fin)
  colnames(fin) <- c("node_ratio", "betweenness", "degree", "closeness", "ivi", "random")
  colnames(fin_p) <- colnames(fin)
  rownames(fin_statistic) <- rownames(g_dat)
  colnames(fin_statistic) <- c("betweenness", "degree", "closeness", "ivi")
  lst <- list(NULL)
  lst$fin <- fin
  lst$fin_p <- fin_p
  lst$fin_statistic <- fin_statistic
  lst$rand <- randomMatrix
  return(lst)
}
