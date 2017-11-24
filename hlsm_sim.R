# This is where hlsm.stan is
setwd("~/Documents/Classes/10-725/hlsm")

########################
# GENERATING HLSM DATA #
########################
# The sparse path connecting the dense ends
w <- cbind(0.5*0:12-3, rep(0,13))
# The two dense barbell ends
x <- mvrnorm(10, mu=c(-3,0), Sigma=array(c(.1,0,0,.1),dim=c(2,2)))
y <- mvrnorm(10, mu=c(3,0), Sigma=array(c(.1,0,0,.1),dim=c(2,2)))
z <- rbind(w,x,y)
plot(z)

rotation.matrix <- function(theta){
    R <- array(c(cos(theta),  sin(theta),
                 -sin(theta), cos(theta)), dim=c(2,2))
    return(function(x){x%*%R})
}
rotate.45 <- rotation.matrix(pi/4)
rotate.90 <- rotation.matrix(pi/2)

plot(rotate.45(z)) # Yup!
plot(rotate.90(z))

N <- nrow(z)
K <- 2
Z <- array(dim=c(N,N,K))

d <- as.matrix(dist(z))
a <- 1

M <- array()


logit <- function(x) 1/(1+exp(-x))

for(i in 1:N){
    for(j in 1:N){
        # This is what a LSM assumes
        X[i,j] <- (logit(a - d[i,j]) > .5)+0
    }
}




# All necessary because of the slightly weird way it's necessary
# to declare data dimensions
sigma  <- 10
sigmat <- array(0, dim=c(K,2,2))
sigmat[,1,1] <- sigma
sigmat[,2,2] <- sigma

hlsm.data <- list(N=N,
                  K=K,
                  edges=edges+0, # without this, stan casts logical => int, which FLATTENS
                  sigma_alpha=10, # NO idea what these should be, using mvlsm's
                  mu_b=c(0,0),
                  sigma_b=sigma*diag(2), # [s   0; 0   s]
                  sigma_z=sigmat)        # [s_k 0; 0 s_k]

# Finding a good starting point
library(MASS) # thought that was automatic
# Sum the edges
summed.edges <- apply(edges, MARGIN=1:2, sum) # THERE; also still symmetric
d <- 1/summed.edges
d[!is.finite(d)] <- 5 # Because. The biggest value in summed.edges = 5.
mds <- isoMDS(as.matrix(d), k=2) # converged! That's good news -- 

# OR EVEN -- use latentnet itself, actually. yeah do that, idiot.
library(latentnet)
# lsm <- ergmm(summed.edges>0 ~ euclidean(d=2))

hlsm.init <- function(){
    # alpha, b, and z; using mds$points
    z <- array(dim=c(N,K,2)) # defined above
    for(k in 1:K)
        z[,k,] <- lsm$mcmc.pmode$Z #mds$points
    list(alpha=rnorm(K), z=z, b=lsm$mcmc.pmode$Z) #mds$points
}

library(rstan)
fit <- stan(file="hlsm.stan",
            model_name="hlsm",
            data=hlsm.data,
            init=hlsm.init,
            verbose=T)
save(fit, file="hlsm_fit.Rdata")

# 1) I want to plot the results (draw network)
# 2) Need to do with TEST data (simple example)
# 3) Classification accuracy sort of thing, pairs
#       accurately represented, or not



