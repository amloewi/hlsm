# --------------------------------
# 36-725 - Fall 2017
# Convex Optmization
# Final Project
# Hierarchical Latent Space Model
# --------------------------------
#

sigmoid <- function(x) {
  # computes sigmoid function
  #
  return(1 / (1 + np.exp(x)))
}

sigmoid.hlsm <- function(b, epsilon, y) {
  # Computes sigmoid function for the HLSM model
  #
  # Args:
  #     b: numeric vector of dimension (1 vs n)
  #     epsilon: numeric array of dimension (n vs k)
  #     y: numeric array of dimension (n vs n vs k)
  # x <- 
}

compute.eta <- function(b, epsilon) {
  # computes the eta function (argument of the sigmoid)
  #
  # Args:
  #     b:  numeric vector of dimension (2 vs n) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (2 vs n vs k) with deviations in 2-d of 
  #               from b for each node in each level

  # check if dimensions will be propagated correctly
  z <- b + epsilon # z should be 2 vs n vs k
  
  # compute L-2 distance
  d <- dist(t(x), diag = TRUE) # d should be n vs n vs k (symmetric matrix)
  
  
  
  
  
  
  
  
  
}

objective.function <- function(b, epsilon, y, lambda) {
  # Computes objective function (Binomial log MLE with lasso penalty)
  #
  # Args:
  #     b: numeric vector of dimension (1 vs n)
  #     epsilon: numeric array of dimension (n vs k)
  #     y: numeric array of dimension (n vs n vs k)
  
  out <- (1-y) * compute.eta(b, epsilon) - log(1 + exp(compute.eta(b, epsilon)))
  out <- sum(out) + lambda * sum(sqrt(epsilon^2))
  
  return(out)
}
