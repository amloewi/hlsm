# ****************************************
# 36-725 - Fall 2017
# Convex Optmization
# Final Project
# Hierarchical Latent Space Model
# ****************************************
#
#
dummy.y <- function(){
  
 y <- array(0, dim=c(10, 10, 3)) 

 y[ , , 1] <- matrix(c(0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                       0, 0, 0, 1, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 1, 0, 0, 0, 1, 0, 0,
                       0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 1, 1, 0, 0, 0, 0, 0,
                       0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
                       0, 0, 1, 1, 1, 0, 0, 1, 0, 0))
   
  
}


test.case <- function(){
  
  N <- 20
  circle <- seq(0,2*pi,length=N)
  b <- cbind(cos(circle), sin(circle))
  
  stretch <- function(x, m, dim="x"){
    stopifnot(any(c("x", "y") %in% dim))
    M <- matrix(c(1, 0, 0, 1), nrow=2)
    if("x" %in% dim) M[1,1] <- m
    if("y" %in% dim) M[2,2] <- m
    return(x %*% M)
  }
  
  stretch(b, 2, dim="x")
  stretch(b, 2, dim="y")
  
  return(b)
  
  
}

sigmoid <- function(x) {
  # computes sigmoid function
  #
  return(1 / (1 + np.exp(x)))
}

compute.eta <- function(b, epsilon, alpha.value=5) {
  # computes the eta function (argument of the sigmoid)
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  
  n <- dim(epsilon)[1]
  K <- dim(epsilon)[3]
  
  alpha <- rep(alpha.value, K)
  
  # check if dimensions will be propagated correctly
  z <- array(b, dim=c(dim(b), K)) + epsilon

  d <- array(0, dim=c(n*(n-1)/2, K))
  # compute square L-2 distance
  # d[ ,k] is a vector of size n*(n-1)/2
  # it is the vector with the elements of lower half of the distance matrix
  # (ordered columnwise - R standard)
  for (k in 1:K){
    d[ ,k] <- as.vector((dist(z[ , , k])^2))
  }
  
  # return a matrix of size (n/2*(n-1) vs K) with the distances
  # within each level
  return(alpha.value - d)
}

objective.function <- function(b, epsilon, y, lambda) {
  # Computes objective function (Binomial log MLE with lasso penalty)
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k)
  
  # because y[, , k] is a symmetric matrix and only need the lower half
  # we create a new matrix y.2 with size (n/2*(n-1) vs K)
  # where each column k is a vector with the elements of the lower half
  # of the matrix y[, , k]
  y.2 <- apply(y, 3, FUN={function(x){as.vector(x[lower.tri(x)])}}) 
  
  out <- (1-y.2) * compute.eta(b, epsilon) - log(1 + exp(compute.eta(b, epsilon)))
  out <- sum(out) + lambda * sum(sqrt(epsilon^2))
  
  return(out)
}

proximal_op <- function(b, epsilon, t, lambda){ 
  # Computes proximal operator
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     t: is a scalar
  #     lambda: is a scalar

  dim.epsilon <- dim(epsilon)

  # out is the resulting array (n vs 2 vs k+1)
  # the first of the k+1 matrices is the b position
  out <- array(0, dim=c(dim(epsilon)[1], dim(epsilon)[2], dim(epsilon)[3]+1))
  out[ , , 1] <- b
  
  for (k in dim.epsilon[3]){
    norm_value <- sqrt(sum(epsilon[ , , k]^2))
    
    if(norm_value >= t*lambda){
      alpha <- ((norm_value - t*lambda)/norm_value)*epsilon[ , , k]
    }
    
    out[ , , k] <- alpha
  }
  
  return(out)
}

gradient_g <- function(b, epsilon, y){
  # Compute the gradient of the MLE function wrt to b and epsilon
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k)
  
  n <- dim(epsilon)[2]
  K <- dim(epsilon)[3]
  
  z <- array(b, dim=c(dim(b), K)) + epsilon
  
  # computes exp(eta)/(1+exp(eta)) and allocates it to a 
  # n vs n vs K array (needed to be coherent with y)
  exp.eta.aux <- exp(compute.eta(b, epsilon))/(1+compute.eta(b, epsilon))
  exp.eta <- array(0, dim=c(n, n, K))
  for(k in 1:K) {
    m <- matrix(0,n,n)
    m[lower.tri(m)] <- exp.eta.aux[ , k]
    exp.eta[ , , k] <- m
  }
  
  # grad is the resulting gradient array (n vs 2 vs k+1)
  # there are K+1 because of the 1 b + K epsilons
  grad <- array(0, dim=c(dim(epsilon)[1], dim(epsilon)[2], dim(epsilon)[3]+1))
  
  # gradient for b
  for (m in 1:n) {
    aux <- 0
    for (k in 1:K) {
      for (j in 1:(m-1)) {
        aux <- aux - (y[m, j, k] - 1 + exp.eta[m, j, k])*(z[m, ,k] - z[j, , k]) 
      }
      
      for (i in (m+1):n) {
        aux <- aux + (y[i, m, k] - 1 + exp.eta[i, m, k])*(z[i, ,k] - z[m, , k])
      }
    }
    grad[m, , 1] <- 2*aux
  }
  
  # gradient for epsilon_k
  for (k in 1:K) {
    aux <- 0
    for (m in 1:n) {
      for (j in 1:(m-1)) {
        aux <- aux - (y[m, j, k] - 1 + exp.eta[m, j, k])*(z[m, ,k] - z[j, , k]) 
      }
      
      for (i in (m+1):n) {
        aux <- aux + (y[i, m, k] - 1 + exp.eta[i, m, k])*(z[i, ,k] - z[m, , k])
      }
    }
    grad[m, , (k+1)] <- 2*aux
  }
  
  # array (n vs 2 vs k+1)
  return(grad)
}

generalized_grad <- function(b, epsilon, y, t, lambda){
  # Computes the generalized gradient for the proximal gradient method
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k)
  #     t: scalar

  n <- dim(epsilon)[1]
  K <- dim(epsilon)[3]
  
  beta <- array(0, dim=c(n, dim(epsilon)[2], K+1))
  beta[ , , 1] <- b
  beta[ , , 2:(K+1)] <- epsilon
  
  grad_g = gradient_g(b, epsilon, y)
  
  beta1 = beta - t*grad_g
  
  gen_grad = (1/t) * (beta - proximal_op(b, epsilon, t, lambda))
  
  return(gen_grad)
}
  
update_beta <- function(b, epsilon, y, t, lambda){
  # Updates the value of the parameter array
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k)
  #     t: scalar
  #     lambda: scalar

  n <- dim(epsilon)[1]
  K <- dim(epsilon)[3]
  
  beta <- array(0, dim=c(n, dim(epsilon)[2], K+1))
  beta[ , , 1] <- b
  beta[ , , 2:(K+1)] <- epsilon
  
  beta_new = beta - t * generalized_grad(b, epsilon, y, t, lambda)
  
  return(beta_new)
}


run.optimization <- function() {
  
  n.steps <- 1000
  
  # read y
  
  n <- dim(y)[1]
  K <- dim(y)[3]
  
  b <- array(0, dim=c(n, 2))
  epsilon <- array(0, dim=c(n, 2, K))
  
  obj.value <- rep(0, n.steps)
  
  for (s in 1:n.steps) {
    # evaluate objective function
    obj.value[s] <- objective.function(b, epsilon, y, lambda)
    cat(sprintf("Step: %4d ; Obj. Value: %15.6f \n", s, obj.value[s]))
    # update parameters
    param <- update_beta(b, epsilon, y, t, lambda)
    b <- param[ , ,1]
    epsilon <- param[ , ,2:(K+1)]
  }
  
  plot(1:n.steps, obj.value, type='l', xlab='step', 
       ylab='Objective Function')
}