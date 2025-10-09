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