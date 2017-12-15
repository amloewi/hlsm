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
  yy <- matrix(0, nrow=10, ncol=10)
  
  aa <- sample(c(0,1), 45, replace=TRUE)
  yy[lower.tri(yy)] <- aa
  yy <- yy + t(yy)
  
  y[ , , 1] <- yy
  
  # randomly delete 4 edges
  aa.2 <- aa
  aa.2[sample(which(aa == 1), 4)] <- 0
  # randomly create 4 more edges
  aa.2[sample(which(aa == 0), 4)] <- 1

  yy <- matrix(0, nrow=10, ncol=10)
  yy[lower.tri(yy)] <- aa.2
  yy <- yy + t(yy)
  
  y[ , , 2] <- yy
  
  # randomly delete 4 edges
  aa.3 <- aa
  aa.3[sample(which(aa == 1), 4)] <- 0
  # randomly create 4 more edges
  aa.3[sample(which(aa == 0), 4)] <- 1
  
  yy <- matrix(0, nrow=10, ncol=10)
  yy[lower.tri(yy)] <- aa.3
  yy <- yy + t(yy)
  
  y[ , , 3] <- yy
  
  g1 <- graph.adjacency(y[ , , 1], mode='undirected')
  plot.igraph(g1, layout=layout.circle)

  g2 <- graph.adjacency(y[ , , 2], mode='undirected')
  plot.igraph(g2, layout=layout.circle)

  g3 <- graph.adjacency(y[ , , 3], mode='undirected')
  plot.igraph(g3, layout=layout.circle)
  
  return(y)
}

sigmoid <- function(x) {
  # computes sigmoid function
  #
  return(1 / (1 + np.exp(x)))
}

compute.eta <- function(b, epsilon, alpha.value=0.3) {
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

log.lik <- function(b, epsilon, y) {
  # because y[, , k] is a symmetric matrix and only need the lower half
  # we create a new matrix y.2 with size (n/2*(n-1) vs K)
  # where each column k is a vector with the elements of the lower half
  # of the matrix y[, , k]
  y.2 <- apply(y, 3, FUN={function(x){as.vector(x[lower.tri(x)])}}) 
  
  out <- (1-y.2) * compute.eta(b, epsilon) + log(1 + exp(-compute.eta(b, epsilon)))
  out <- sum(out)
  
  return(out)
}

objective.function <- function(b, epsilon, y, lambda, group.level=FALSE) {
  # Computes objective function (Binomial log MLE with lasso penalty)
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k)
  #     group.level: (logical) should we group penalties on position by level?
  
  # because y[, , k] is a symmetric matrix and only need the lower half
  # we create a new matrix y.2 with size (n/2*(n-1) vs K)
  # where each column k is a vector with the elements of the lower half
  # of the matrix y[, , k]
  y.2 <- apply(y, 3, FUN={function(x){as.vector(x[lower.tri(x)])}}) 
  
  out <- (1-y.2) * compute.eta(b, epsilon) + log(1 + exp(-compute.eta(b, epsilon)))
  
  if (group.level) {
    # grouped lasso (like hw3). All groups have same weight
    out <- sum(out) + lambda * sum(sqrt(epsilon^2))
  } else {
    # pointwise lasso. Use L-1 norm
    out <- sum(out) + lambda * sum(abs(epsilon))
  }
  
  return(out)
}

proximal_op <- function(b, epsilon, t, lambda, group.level=FALSE){ 
  # Computes proximal operator for grouped and non grouped cases
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with positions in 2-d of 
  #         each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with deviations in 
  #               2-d of from b for each node in each level
  #     t: is a scalar
  #     lambda: is a scalar
  #     group.level: (logical) should we group penalties on position by level?

  dim.epsilon <- dim(epsilon)

  # out is the resulting array (n vs 2 vs k+1)
  # the first of the k+1 matrices is the b position
  out <- array(0, dim=c(dim(epsilon)[1], dim(epsilon)[2], dim(epsilon)[3]+1))
  out[ , , 1] <- b
  
  if (group.level) {
    # penalties on deviations grouped by level. Use L-2 Norm (like HW3)
    for (k in 1:dim.epsilon[3]){
      norm_value <- sqrt(sum(epsilon[, , k]^2))
      alpha <- 0
      if(norm_value >= t*lambda){
        alpha <- ((norm_value - t*lambda)/norm_value)*epsilon[, , k]
      }
      out[, , k+1] <- alpha
      
      # out[, , k+1] <- sign(epsilon[, , k])*pmax(abs(epsilon[, , k]) - t*lambda, 0)
    }
  } else {
    # penalties on deviations are NOT grouped by level. Use L-1 Norm
    for (k in 1:dim.epsilon[3]){
      for (i in 1:dim.epsilon[1]) {
        # norm_value <- sqrt(sum(epsilon[i, , k]^2))
        # alpha <- 0
        # if(norm_value >= t*lambda){
        #  alpha <- ((norm_value - t*lambda)/norm_value)*epsilon[i, , k]
        # }
        # out[i, , k+1] <- alpha
        
        out[i, , k+1] <- sign(epsilon[i, , k])*pmax(abs(epsilon[i, , k]) - t*lambda, 0)
      }
    }
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
  
  n <- dim(epsilon)[1]
  K <- dim(epsilon)[3]
  
  z <- array(b, dim=c(dim(b), K)) + epsilon
  
  # computes exp(eta)/(1+exp(eta)) and allocates it to a 
  # n vs n vs K array (needed to be coherent with y)
  exp.eta.aux <- exp(-compute.eta(b, epsilon))/(1+exp(-compute.eta(b, epsilon)))
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
    aux <- c(0, 0)
    for (k in 1:K) {
      if(m > 1){
        for (j in 1:(m-1)) {
          aux <- aux - (y[m, j, k] - 1 + exp.eta[m, j, k])*(z[m, ,k] - z[j, , k]) 
        }
      }
      if(m<n){
        for (i in (m+1):n) {
          aux <- aux + (y[i, m, k] - 1 + exp.eta[i, m, k])*(z[i, ,k] - z[m, , k])
        }
      }
    }
    # the negative sign here comes from the fact that it is the
    # negative log likelihood
    grad[m, , 1] <- -2*aux
  }
  
  # gradient for epsilon_k
  for (k in 1:K) {
    for (m in 1:n) {
      aux <- c(0, 0)
      if(m>1){
        for (j in 1:(m-1)) {
          aux <- aux - (y[m, j, k] - 1 + exp.eta[m, j, k])*(z[m, ,k] - z[j, , k]) 
        }
      }
      
      if(m<n){
        for (i in (m+1):n) {
          aux <- aux + (y[i, m, k] - 1 + exp.eta[i, m, k])*(z[i, ,k] - z[m, , k])
        }
      }
      # the negative sign here comes from the fact that it is the
      # negative log likelihood
      grad[m, , (k+1)] <- -2*aux
    }
  }
  
  # array (n vs 2 vs k+1)
  return(grad)
}

update_beta <- function(b, epsilon, y, t, lambda, group.level=FALSE){
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
  
  grad_g <- gradient_g(b, epsilon, y)
  
  beta.aux <- beta - t * grad_g
  b.aux <- beta.aux[ , , 1]
  epsilon.aux <- beta.aux[ , , 2:(K+1)]
  
  beta.new <- proximal_op(b.aux, epsilon.aux, t, lambda, group.level=group.level)
  
  return(beta.new)
}

run.optimization <- function(b, epsilon, y, lambda=1, t=1e-3, 
                             group.level=FALSE, n.steps=3000) {
  # Updates the value of the parameter array
  #
  # Args:
  #     b:  numeric vector of dimension (n vs 2) with initial positions  
  #         in 2-d of each node
  #     epsilon:  numeric array of dimension (n vs 2 vs k) with initial
  #               deviations in 2-d of from b for each node in each level
  #     y: numeric array of dimension (n vs n vs k) with adjency matrix in each
  #        level  
  #     t: (scalar) step size
  #     lambda: (scalar) weight of penalty term in optimization function
  #     group.level: (logical) should we group penalties on position by level?
  #     n.steps: (integer) number of iterations in the optimization loop
  #               obs: we could change this to a stopping criteria, of the kind
  #                    | f^(k) - f^(k-1) | < tolerance
  
  n <- dim(y)[1]
  K <- dim(y)[3]
  
  obj.value <- rep(0, n.steps)
  
  for (s in 1:n.steps) {
    # evaluate objective function
    obj.value[s] <- objective.function(b, epsilon, y, lambda, group.level=group.level)
    cat(sprintf("Step: %4d ; Obj. Value: %15.6f \n", s, obj.value[s]))
    # update parameters
    param <- update_beta(b, epsilon, y, t, lambda, group.level=group.level)
    b <- param[ , ,1]
    epsilon <- param[ , ,2:(K+1)]
  }
  
  # plot objective function in standard output
  plot(1:n.steps, obj.value, type='l', xlab='step', 
       ylab='Objective Function')
  
  # return list with results
  return(list(b=b, epsilon=epsilon, obj.value=obj.value))
}

compute.z <- function(b, epsilon) {
  K <- dim(epsilon)[3]
  z <- array(b, dim=c(dim(b), K)) + epsilon
  return(z)
}

run.opt.wrapper <- function(){
  # uses alex initial values to call optimization
  
  source('../hlsm_sim_francisco.R')
  
  list.results.total <- list()
  for (case in c('individual', 'group')) {
    b <- ovals.init$b
    n <- dim(b)[1]
    
    # ovals.init$z is (n,K,2)
    z <- ovals.init$z
    z <- aperm(z, c(1, 3, 2)) # permute to (n, 2, K)
    K <- dim(z)[3]
    
    epsilon <-  z - array(b, dim=c(dim(b), K))
    y <- two.ovals
    
    gf <- case == 'group'
    list.results <- list()
    lambda.values <- c(0.1, 1, 10, 100)
    for (i in 1:length(lambda.values)){
      result.opt <- run.optimization(b, epsilon, y, lambda = lambda.values[i], 
                                     group.level = gf)
      list.results[[i]] <- result.opt
    }
    
    figure.width <- 0.8 * 8.5 # inches
    
    library(ggplot2)
    library(plyr)
    library(dplyr)
    library(tidyr)
    
    n.steps <- length(list.results[[1]]$obj.value)
    
    df.obj.values <- data.frame(matrix(NA,nrow=n.steps, 
                                       ncol=length(lambda.values)))
    names(df.obj.values) <- lambda.values
    for (i in 1:length(lambda.values)) {
      df.obj.values[,i] <- list.results[[i]]$obj.value
    }
    df.obj.values$step <- 1:nrow(df.obj.values)
    
    df.obj.values <- gather(df.obj.values, lambda, value, `0.1`:`100`, factor_key=TRUE)
    
    g1 <- ggplot() + geom_line(data=df.obj.values, 
                               aes(x=step, y=value, linetype=lambda))
    
    g1 <- g1 + theme_bw() + scale_y_continuous(trans='log') + 
      theme(axis.text.y = element_blank(), legend.position=c(0.95,0.95), 
            legend.justification=c(1,1)) + ylab('Log Obj Function') + 
      xlab('iteration')
    
    nf <- paste0('plot_objective_', case, '.png')
    png(nf, res=300, width=figure.width, height=3/5*figure.width, unit="in")
    print(g1)
    dev.off()
    
    nf <- paste0('plot_positions_', case, '.png')
    png(nf, res=300, width=figure.width, height=figure.width, unit="in")
    par(mfrow=c(2, 2), mar=c(0, 0, 1, 0)+0.1)
    plot.positions(z, b, y)
    title(main='Initial Positions')
    for (i in 1:(length(lambda.values)-1)) {
      b.new <- list.results[[i]]$b
      z.new <- array(list.results[[i]]$b, dim=c(dim(list.results[[i]]$b), K)) + 
        list.results[[i]]$epsilon
      plot.positions(z.new, b.new, y)
      title(main=bquote(lambda == .(lambda.values[i])))
    }
    dev.off()
    list.results.total[[case]] <- list.results
  }
  saveRDS(list.results.total, file='complete_results.rds')
}




# ************************************************************
# ***** TRASH ******---- 
# while (TRUE) {
#   gen_grad <- generalized_grad(b, epsilon, y, t.step, lambda)
#   beta.new <- beta - t.step * gen_grad
#   b.new <- beta.new[ , , 1]
#   epsilon.new <- beta.new[ , , 2:(K+1)]
#   
#   if(log.lik(b.new, epsilon.new, y) > log.lik(b, epsilon, y)
#      - t.step * sum(grad_g * gen_grad) + t.step/2 * sum(gen_grad^2)){
#     t.step <- t.step * beta.step
#   } else {
#     break
#   }
#   
#   if(ii > 100) {
#     cat('back tracking line search ended without gowing down...\n')
#     break
#   }
#   ii <- ii+1
# }

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
  
  grad_g <- gradient_g(b, epsilon, y)
  
  beta1 <- beta - t*grad_g
  
  b1 <- beta1[ , , 1]
  epsilon1 <- beta1[ , , 2:(K+1)]
  
  gen_grad = (1/t) * (beta - proximal_op(b1, epsilon1, t, lambda))
  
  return(gen_grad)
}


