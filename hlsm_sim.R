# This is where hlsm.stan is
setwd("~/Documents/Research/hlsm")

######################################
# FUNCTIONS FOR GENERATING HLSM DATA #
######################################

# Linearly scales a set of data points by the scalar 'm', along the axis 'dim'(="x"/"y")
stretch <- function(x, m, dim="x"){
    stopifnot(any(c("x", "y") %in% dim))
    M <- matrix(c(1, 0, 0, 1), nrow=2)
    if("x" %in% dim) M[1,1] <- m
    if("y" %in% dim) M[2,2] <- m
    return(x %*% M)
}

# # Looks good
####### FUNCTION TESTS #######
# plot(b, xlim=c(-2,2), ylim=c(-2,2))
# points(stretch(b, 2, dim=c("x", "y")), col="red")
# points(stretch(b, 2, dim="y"), col="blue")
# plot.lsm(1, b, add=T)
# plot.lsm(1, stretch(b, 2), add=T, col="blue")
# plot.lsm(1, stretch(b, 2, dim="y"), add=T, col="red")
# plot.lsm(1, rotate.45(stretch(b, 2, dim="y")), add=T, col="green")



# Closure. Takes an angle (in radians) and returns a function that rotates by that angle.
rotation.matrix <- function(theta){
    R <- array(c(cos(theta),  sin(theta),
                 -sin(theta), cos(theta)), dim=c(2,2))
    return(function(x){x%*%R})
}
rotate.45 <- rotation.matrix(pi/4)
rotate.90 <- rotation.matrix(pi/2)

####### FUNCTION TESTS #######
# plot(rotate.45(z)) # Yup!
# plot(z, ylim=c(-5,5))
# points(rotate.90(z))


logit <- function(x) 1/(1+exp(-x))

# Takes the parameters of a latent space model and produces an adjacency matrix
positions.to.matrix <- function(a, z){

    N <- nrow(z)
    d <- as.matrix(dist(z))^2 # The SQUARED euclidean distances, remember ... (right?)
    M <- array(dim=c(N,N))
    
    for(i in 1:N){
        for(j in 1:N){
                M[i,j] <- (logit(a - d[i,j]) > .5)+0
        }
    }
    return(M)
}

# 
positions.to.multigraph <- function(Z, alpha=1){
    N <- nrow(Z[[1]])
    K <- length(Z)
    edges <- array(dim=c(N,N,K))
    for(i in 1:K){
        edges[,,i] <- positions.to.matrix(Z[[i]], a=alpha)
    }
    return(edges)
}

# (Functions for ...)
#################################
# Finding a good starting point #
#################################

# hlsm.init <- function(){
#     list(alpha=c(1,1),
#          z=abind(b,b,along=3), # initializes the ovals as circles, also
#          b=b)
# }

find.init <- function(m){
    N <- dim(m)[1]
    K <- dim(m)[3]
    avg <- (apply(m, MARGIN=1:2, sum) > 0) + 0    # init failed with K=3; round=>0?
    lsm <- ergmm(avg ~ euclidean(d=2)) # Squared?
    z.hat <- lsm$mcmc.pmode$Z
    init <- function(){
        z <- array(dim=c(N,K,2))
        for(k in 1:K)
            z[,k,] <- z.hat
        return(list(alpha=rep(1,K), z=z, b=z.hat))
    }
    return(init) # needs a FUNCTION, remember -- why, who knows.
}

# (Functions for ...)
################
# RUNNING STAN #
################

library(latentnet)


# Designed to go from parameters to model in one line
one.shot <- function(edges, sigma, file, model_name, iter){
    
    N <- dim(edges)[1]
    K <- dim(edges)[3]
    
    if(file=="hlsm.stan"){
        sigmat <- array(0, dim=c(K,2,2))
        sigmat[,1,1] <- sigma
        sigmat[,2,2] <- sigma
    } else {
        sigmat <- array(sigma, dim=c(K,2))
    }
    
    hlsm.data <- list(N=N,
                      K=K,
                      edges=edges,     # make sure these are ints  (not bool) else FLATTENED
                      sigma_alpha=10,   # NO idea what these should be, using mvlsm's
                      mu_b=c(0,0),
                      sigma_b=sigma*diag(2), # [s   0; 0   s] (one matrix)
                      sigma_z=sigmat)        # [s_k 0; 0 s_k] (one matrix per layer)
    
    hlsm.init <- find.init(edges)
    fit <- stan(file=file,
                model_name=model_name,
                iter=iter,
                data=hlsm.data,
                init=hlsm.init,
                verbose=T,
                control=list(max_treedepth = 15))
    theta <- unpack.hlsm.fit(fit, N, K)
    return(list(fit=fit, theta=theta, init=hlsm.init()))
}



# (Functions to ...)
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

# zhat <- unpack.hlsm.fit(fit, N, K)
# plot.multigraph(zhat$z, alpha=.02)


# Takes the parameters for an lsm, lays out a network as specified
# WOULD LIKE TO HAVE A FEATURE WHERE -- black lines for actual, red for mistaken edges
plot.lsm <- function(alpha, z, add=F, col="black", xlim=c(-2,2), ylim=c(-2,2)){
    edges <- positions.to.matrix(alpha, z)==1 # => logical
    if(add){
        points(z, col=col)
    } else {
        plot(z, col=col, xlim=xlim, ylim=ylim)
    }
    N <- nrow(z)
    for(i in 1:N){
        for(j in 1:N){
            if(edges[i,j]){
                lines(z[c(i,j),1], z[c(i,j),2], col=col)
            }
        }
    }
}

plot.multigraph <- function(Z, alpha=1){
    plot.lsm(alpha, Z[,1,], add=F, col=1) # WATCH THE INDEXING HERE -- DID I CHANGE THIS?
    for(i in 2:dim(Z)[3]){
        plot.lsm(1, Z[,2,], add=T, col=i)
    }
}

plot.model <- function(m, alpha=NULL, xlim=NULL, ylim=NULL){
    z <- m$theta$z
    b <- m$theta$b
    if(is.null(alpha))
        alpha <- m$theta$alpha
    # First, make sure the plot is big enough
    if(is.null(xlim))
        xlim <- range(c(z[,,1], b[,1]))
    if(is.null(ylim))
        ylim <- range(c(z[,,2], b[,2]))
    
    plot.lsm(alpha[1], z[,1,], add=F, col=1, xlim=xlim, ylim=ylim)
    for(i in 2:dim(z)[2]){
        plot.lsm(alpha[i], z[,i,], add=T, col=i)
    }
}

####### FUNCTION TESTS #######
# plot.multigraph(zhat$z, alpha=.02)
# 
# plot.lsm(1, b)
# plot(b, xlim=c(-.5,.5), ylim=c(-.5,.5))
# plot.lsm(.02, zhat$b,      add=T)
# plot.lsm(.02, zhat$z[,1,], add=T, col="red") # .01 => disconnected
# plot.lsm(.02, zhat$z[,2,], add=T, col="blue")




#######################################################
# SETTING PARAMETERS CREATING DATA AND RUNNING MODELS #
#######################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # => 4 chains at once. Issues?


############
# THE DATA #
############

# library(igraph)
# out <- positions.to.matrix(stretch(b, 3), a=1) # moderate clustering at 2, serious at 3
# plot(graph_from_adjacency_matrix(out)) # sure, not bad

two.circles <- positions.to.multigraph(list(stretch(b, 2, dim=c("x", "y")),
                                            b))

# For the L1 test -- lambda?
two.overlapping.circles   <- positions.to.multigraph(list(b, b))
three.overlapping.circles <- positions.to.multigraph(list(b, b, b))


# # # # # # # # # # # # # # # 
##############################
# RUNNING AND TESTING GRAPHS #
##############################

N <- 20
circle <- seq(0,2*pi,length=N)
b <- cbind(cos(circle), sin(circle))
niter <- 1e3
sigma <- 10

# TWO ORTHOGONAL OVALS
two.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="x"),
                                          stretch(b, 2, dim="y")))
two.ovals.model <- one.shot(two.ovals, sigma, "hlsm.stan", "two_ovals", niter)
plot.multigraph(two.ovals.model$theta$z)
plot.model(two.ovals.model)


# THREE ~ORTHOGONAL OVALS
three.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="x"),
                                            stretch(b, 2, dim="y"),
                                            rotate.45(stretch(b, 2, dim="y"))))
three.ovals.model <- one.shot(three.ovals, sigma, "hlsm.stan", "three_ovals", 2e4)
plot.model(three.ovals.model, alpha=c(.10, .10, .10), xlim=c(-1,1), ylim=c(-1,1)) # oh VERY interesting -- what the hell?


# TWO CONCENTRIC CIRCLES
two.circles.model <- one.shot(two.circles, sigma, "hlsm.stan", "two_circles", niter)

two.ovlp.circles.model   <- one.shot(two.overlapping.circles,  sigvec, "hlsm_lasso.stan", "two_ovlp_circles",    1e3)
#three.ovlp.circles.model <- one.shot(three.ovlp.circles, sigvec, "hlsm_lasso.stan", "three_ovlp__circles", 2e4)









# What's ONE of these -- 
hlsm.data[["edges"]]   <- two.circles
hlsm.data[["sigma_z"]] <- sigmat
hlsm.init <- find.init(two.circles)
two.circles.fit <- stan(file="hlsm.stan",
                        model_name="hlsm_sim_two_circles",
                        iter=1e4,
                        data=hlsm.data, init=hlsm.init,
                        verbose=T, control=list(max_treedepth = 15))
# BEWARE: These are only defined BELOW. Will need to shuffle everything.
zhat <- unpack.hlsm.fit(two.circles.fit, N, K)
plot.multigraph(zhat$z, alpha=5)







#########################################
#########################################
# THE DUMP -- old tests, bad ideas, etc #
#########################################
#########################################

# BARBELLS AND ROTATIONS -- maybe not very clear thinking
# The sparse path connecting the dense ends
w <- cbind(0.5*0:12-3, rep(0,13))
# The two dense barbell ends
x <- mvrnorm(10, mu=c(-3,0), Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
y <- mvrnorm(10, mu=c(3,0),  Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
z <- rbind(w,x,y)
plot(z)


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

hlsm.init <- function(){
    # alpha, b, and z; using mds$points
    z <- array(dim=c(N,K,2)) # defined above
    for(k in 1:K)
        z[,k,] <- lsm$mcmc.pmode$Z #mds$points
    list(alpha=rnorm(K), z=z, b=lsm$mcmc.pmode$Z) #mds$points
}



# ALL THIS WAS A GOOD FIRST TEST, BUT I DON'T NEED IT NOW --
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



fit.20k <- stan(file="hlsm.stan",
             model_name="hlsm_sim",
             data=hlsm.data,
             init=hlsm.init,
             iter=2e4,
             verbose=T,
             control=list(max_treedepth = 15)) # from 10, as suggested
# 2k was lousy, 5k better-looking (complete eyeball of the traceplot) 20k -- LOOKS ok
# fit.2k <- fit
# fit <- fit.20k
# save(fit, file="hlsm_sim_fit_20k.Rdata")
# zhat.table <- summary(fit)$summary[,"mean", drop=F]

# And finally, plot to check
ghat.1 <- positions.to.matrix(.01, zhat$z[,1,])
ghat.2 <- positions.to.matrix(.01, zhat$z[,2,])
ghat.b <- positions.to.matrix(.01, zhat$b)

# plot(graph_from_adjacency_matrix(out)) # it's -- a big fat ball. Okay both are.
# WAIT -- interesting -- if I just MAKE alpha=1, I get ~ the right answer
# (if, remember, I also initialized it that way. Wait did I?)
plot(graph_from_adjacency_matrix(positions.to.matrix(1, b)))

# ~Classification accuracies
# What happened to 'edges'? It's not a table anymore
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


