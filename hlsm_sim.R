# This is where hlsm.stan is
setwd("~/Documents/Research/hlsm")

########################
# GENERATING HLSM DATA #
########################
N <- 20
K <- 2
circle <- seq(0,2*pi,length=N)
b <- cbind(cos(circle), sin(circle))

# Linearly scales a set of data points by the scalar 'm', along the axis 'dim'(="x"/"y")
stretch <- function(x, m, dim="x"){
    stopifnot(dim=="x" | dim=="y")
    if(dim=="x") M <- matrix(c(m, 0, 0, 1), nrow=2)
    if(dim=="y") M <- matrix(c(1, 0, 0, m), nrow=2)
    return(x %*% M)
}

# Looks good
plot(b, xlim=c(-2,2), ylim=c(-2,2))
points(stretch(b, 2))
points(stretch(b, 2, dim="y"))


logit <- function(x) 1/(1+exp(-x))

# Takes the parameters of a latent space model and produces an adjacency matrix
positions.to.matrix <- function(a, z){

    N <- nrow(z)
    d <- as.matrix(dist(z))^2 # The SQUARED euclidean distances, remember
    M <- array(dim=c(N,N))
    
    for(i in 1:N){
        for(j in 1:N){
                M[i,j] <- (logit(a - d[i,j]) > .5)+0
        }
    }
    return(M)
}


library(igraph)
out <- positions.to.matrix(stretch(b, 3), a=1) # moderate clustering at 2, serious at 3
plot(graph_from_adjacency_matrix(out)) # sure, not bad

############
# THE DATA #
############
# Two ovals of data, at 90 degrees to one another. Hope to infer a circle of base positions
edges <- array(dim=c(N,N,K))

# How do you stack 2d structures? ALONG=3 (below)
# https://stackoverflow.com/questions/14857358/stacking-array-slices-rbind-for-2-dimensions
library(abind)
# Again, this isn't ... fuckin, exactly like the DGP, is it. Fuckin' ... 
positions <- abind(stretch(b, 3, dim="x"),
                   stretch(b, 3, dim="y"), along = 3)
for(i in 1:2){
    edges[,,i] <- positions.to.matrix(positions[,,i], a=1)
}



################################
# SETTING PRIORS AND STAN DATA #
################################

# All necessary because of the slightly weird way it's necessary
# to declare data dimensions
sigma  <- 10
sigmat <- array(0, dim=c(K,2,2))
sigmat[,1,1] <- sigma
sigmat[,2,2] <- sigma

hlsm.data <- list(N=N,
                  K=K,
                  edges=edges,    # make sure these are ints  (not bool) else FLATTENED
                  sigma_alpha=10, # NO idea what these should be, using mvlsm's
                  mu_b=c(0,0),
                  sigma_b=sigma*diag(2), # [s   0; 0   s] (one matrix)
                  sigma_z=sigmat)        # [s_k 0; 0 s_k] (one matrix per layer)

#################################
# Finding a good starting point #
#################################
# library(latentnet) # USE LATENT NET, LAYER BY LAYER
# lsm <- ergmm(summed.edges>0 ~ euclidean(d=2))

hlsm.init <- function(){
    list(alpha=c(1,1),
         z=abind(b,b,along=3), # initializes the ovals as circles, also
         b=b)
}

################
# RUNNING STAN #
################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # => 4 chains at once. Issues?

fit <- stan(file="hlsm.stan",
            model_name="hlsm_sim",
            data=hlsm.data,
            init=hlsm.init,
            #iter=5e3,
            verbose=T)
save(fit, file="hlsm_sim_fit.Rdata")
zhat.table <- summary(fit)$summary[,"mean", drop=F]

################################
# EXAMINE THE ESTIMATED VALUES #
################################

unpack.hlsm.fit <- function(fit, N, K){
    s <- summary(fit)$summary[,"mean"]
    params <- t(array(s, dim=c(K, 1+N+N*K))) # => one k-dim'nl param per line, CUTS lp__
    return(list(alpha=params[1,],
                b=params[2:(N+1),],
                # => z[N,K,2], i.e. nodes, layers, x/y coords.
                z=aperm(array(params[(N+2):(1+N+N*2),], # This is UGLY, but works
                              dim=c(K,N,2)),
                        c(2,1,3)))) 
}

zhat <- unpack.hlsm.fit(fit, N, K)

# And finally, plot to check
k <- 1
ghat.1 <- positions.to.matrix(1, zhat$z[,1,])
ghat.2 <- positions.to.matrix(1, zhat$z[,2,])
ghat.b <- positions.to.matrix(1, zhat$b)

plot(graph_from_adjacency_matrix(out)) # it's -- a big fat ball. Okay both are.
# WAIT -- interesting -- if I just MAKE alpha=1, I get ~ the right answer
# (if, remember, I also initialized it that way. Wait did I?)
plot(graph_from_adjacency_matrix(positions.to.matrix(1, b)))

# ~Classification accuracies
table(as.vector(ghat.1), as.vector(edges[,,1])) # NOT TERRIBLE!
table(as.vector(ghat.2), as.vector(edges[,,2])) # NOT TERRIBLE! ALSO!
table(as.vector(ghat.b), as.vector(positions.to.matrix(b, a=1))) # Seriously not bad
# Now try tuning with alpha, maybe? And why was it so high?
# Also, try with initialization values that aren't cheating.


out <- array(dim=20)
vals <- seq(0, 5, length=20)
for(i in 1:20){
    cm <- table(as.vector(ghat.b), as.vector(positions.to.matrix(b, a=vals[i]))) 
    out[i] <- sum(diag(cm))/sum(cm)
}
plot(out, xaxt="n")
axis(1, at=1:20, labels=round(vals,2)) # Okay so yes, peaks around 1; why estimate=20?







#########################################
# THE DUMP -- ideas that didn't make it #
#########################################

# BARBELLS AND ROTATIONS -- maybe not very clear thinking
# The sparse path connecting the dense ends
w <- cbind(0.5*0:12-3, rep(0,13))
# The two dense barbell ends
x <- mvrnorm(10, mu=c(-3,0), Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
y <- mvrnorm(10, mu=c(3,0),  Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
z <- rbind(w,x,y)
plot(z)

# Closure. Takes an angle (in radians) and returns a function that rotates by that angle.
rotation.matrix <- function(theta){
    R <- array(c(cos(theta),  sin(theta),
                 -sin(theta), cos(theta)), dim=c(2,2))
    return(function(x){x%*%R})
}
rotate.45 <- rotation.matrix(pi/4)
rotate.90 <- rotation.matrix(pi/2)

plot(rotate.45(z)) # Yup!
plot(z, ylim=c(-5,5))
points(rotate.90(z))

#####################################
# Finding a good starting point 0.1 #
#####################################
# Better idea was to use latentnet itself -- MDS necessary, if you don't HAVE that.
library(MASS) # thought that was automatic (loads isoMDS)
# Sum the edges
summed.edges <- apply(edges, MARGIN=1:2, sum) # THERE; also still symmetric
d <- 1/summed.edges
d[!is.finite(d)] <- 5 # Because. The biggest value in summed.edges = 5.
mds <- isoMDS(as.matrix(d), k=2) # converged! That's good news -- 

# hlsm.init <- function(){
#     # alpha, b, and z; using mds$points
#     z <- array(dim=c(N,K,2)) # defined above
#     for(k in 1:K)
#         z[,k,] <- lsm$mcmc.pmode$Z #mds$points
#     list(alpha=rnorm(K), z=z, b=lsm$mcmc.pmode$Z) #mds$points
# }


