stability_selection_split_data <- function(X, dim = c('rows', 'columns')){

  if (dim == 'rows'){
    n <- nrow(X)
    n1 <- ceiling(n/2)
    # n2 <- floor(nrow(X)/2)
  } else if (dim == 'columns'){
    n <- ncol(X)
    n1 <- ceiling(n/2)
    # n2 <- floor(ncol(X)/2)
  } else {
    stop('Wrong input for dim')
  }

  subset1 <- sample(n, size = n1, replace = FALSE)
  if (dim == 'columns'){
    X1 <- X[,subset1]
    X2 <- X[,-(subset1)]
  } else { # dim = 'rows'
    X1 <- X[subset1,]
    X2 <- X[-(subset1),]
  }
  return(list(X1, X2, subset1, setdiff(c(1:n), subset1)))
}

compute_cosine_sim_matrix <- function(L1, L2){
  norms1 <- apply(L1, 2, function(x){sqrt(sum(x^2))})
  norms1[norms1 == 0] <- Inf

  norms2 <- apply(L2, 2, function(x){sqrt(sum(x^2))})
  norms2[norms2 == 0] <- Inf

  L1_normalized <- t(t(L1)/norms1)
  L2_normalized <- t(t(L2)/norms2)

  #compute matrix of cosine similarities
  cosine_sim_matrix <- crossprod(L1_normalized, L2_normalized)

  return(cosine_sim_matrix)
}

stability_selection_post_processing <- function(L1, L2, threshold=0.99){
  # use cosine similarity as similarity metric
  # could try correlation but may run into issues if there is a constant baseline

  # compute cosine similarity
  cosine_sim_matrix <- compute_cosine_sim_matrix(L1, L2)

  #compute maximum similarity for each factor
  idx <- which.min(dim(cosine_sim_matrix))

  #pick factors to keep
  max_similarity <- apply(cosine_sim_matrix, idx, max)
  factors.idx.keep <- which(max_similarity > threshold)
  if (length(factors.idx.keep) == 0){
    L_final = matrix(0, nrow = nrow(L1), ncol = 1)
  } else if (idx == 1) {
    L_final <- L1[,factors.idx.keep, drop = FALSE]
  } else { # captures idx == 2
    L_final <- L2[, factors.idx.keep, drop = FALSE]
  }
  return(list(L = L_final, similarity = cosine_sim_matrix))
}

stability_selection_row_split_post_processing <- function(L1, L2, F1, F2,
                                                          threshold=0.99){
  # use cosine similarity as similarity metric
  # could try correlation but may run into issues if there is a constant baseline
  n <- nrow(L1) + nrow(L2)

  # compute cosine similarity
  cosine_sim_matrix <- compute_cosine_sim_matrix(F1, F2)

  #compute maximum similarity for each factor
  idx <- which.min(dim(cosine_sim_matrix))

  #pick factors to keep
  max_similarity <- apply(cosine_sim_matrix, idx, max)
  max_similarity_est2_idx <- apply(cosine_sim_matrix, idx, which.max)
  factors.idx.keep <- which(max_similarity > threshold)

  # option 1: concatenate the L estimates; doesn't feel the best?
  if (length(factors.idx.keep) == 0){
    L_final = matrix(0, nrow = n, ncol = 1)
  } else if (idx == 1) {
    print('idx = 1')
    L_final <- rbind(L1[,factors.idx.keep, drop = FALSE], L2[,max_similarity_est2_idx[factors.idx.keep], drop = FALSE])
  } else { # captures idx == 2
    print('idx = 2')
    L_final <- rbind(L1[,max_similarity_est2_idx[factors.idx.keep]], L2[, factors.idx.keep])
  }

  # option 2: get an average estimate for F and then get an estimate for L using least squares or something
  return(list(L = L_final, similarity = cosine_sim_matrix))
}
