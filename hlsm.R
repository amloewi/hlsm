setwd("~/Documents/Research/gsa/study/data")
# edges <- read.csv("anon.edges")
load("staff.Rdata") # the variable I want is 'networks', although it's still dirty
# This is where hlsm.stan is
setwd("~/Documents/Classes/10-725/hlsm")

responded <- !is.na(networks[[1]][,1]) # missing responders
on.survey <- !is.na(networks[[1]][1,]) # 
ok <- responded & on.survey    
N <- sum(ok)
K <- length(networks)
edges <- array(dim=c(N,N,K))

library(sna) # for 'symmetrize'
for(k in 1:K){
    edges[,,k] <- symmetrize(networks[[k]][ok,ok] > 0, rule="weak")
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



