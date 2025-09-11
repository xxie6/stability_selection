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
  X1 <- X[,subset1]
  X2 <- X[,-(subset1)]
  return(list(X1, X2))
}

stability_selection_post_processing <- function(L1, L2, threshold=0.99){
  # use cosine similarity as similarity metric
  # could try correlation but may run into issues if there is a constant baseline
  norms1 <- apply(L1, 2, function(x){sqrt(sum(x^2))})
  norms1[norms1 == 0] <- Inf
  
  norms2 <- apply(L2, 2, function(x){sqrt(sum(x^2))})
  norms2[norms2 == 0] <- Inf
  
  L1_normalized <- t(t(L1)/norms1)
  L2_normalized <- t(t(L2)/norms2)
  
  #compute matrix of cosine similarities
  cosine_sim_matrix <- abs(crossprod(L1_normalized, L2_normalized))
  
  #compute maximum similarity for each factor
  idx <- which.min(dim(cosine_sim_matrix))
  
  #pick factors to keep
  max_similarity <- apply(cosine_sim_matrix, idx, max)
  factors.idx.keep <- which(max_similarity > threshold)
  if (length(factors.idx.keep) == 0){
    return(list(L = matrix(0, nrow = nrow(L1), ncol = 1), 
                similarity = cosine_sim_matrix))
  }
  # I think this could lead to redundant factors?
  # do I want to switch it to be a 1-1 mapping?
  
  # post-processing step to average both estimates?
  
  if(idx == 1){
    L_final <- L1[,factors.idx.keep]
  } else {
    L_final <- L2[, factors.idx.keep]
  }
  return(list(L = L_final, similarity = cosine_sim_matrix))
}