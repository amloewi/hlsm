# This is where hlsm.stan is
setwd("~/Documents/Research/hlsm")


            ##### FUNCTIONS FOR GENERATING HLSM DATA #####

# Linearly scales a set of data points by the scalar 'm', along the axis 'dim'(="x"/"y")
stretch <- function(x, m, dim="x"){
    stopifnot(any(c("x", "y") %in% dim))
    M <- matrix(c(1, 0, 0, 1), nrow=2)
    if("x" %in% dim) M[1,1] <- m
    if("y" %in% dim) M[2,2] <- m
    return(x %*% M)
}

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
rotate <- function(x, theta){
    f <- rotation.matrix(theta)
    return(f(x))
}

# rotate.45 <- rotation.matrix(pi/4)
# rotate.60 <- rotation.matrix(pi/3)
# rotate.90 <- rotation.matrix(pi/2)

####### FUNCTION TESTS #######
# plot(rotate.45(z)) # Yup!
# plot(z, ylim=c(-5,5))
# points(rotate.90(z))


invlogit <- function(x) 1/(1+exp(-x))

# Takes the parameters of a latent space model and produces an adjacency matrix
positions.to.matrix <- function(a, z){

    N <- nrow(z)
    d <- as.matrix(dist(z))^2 # The SQUARED euclidean distances, remember ... (right?)
    M <- array(dim=c(N,N))
    #if(length(a)!=)
    for(i in 1:N){
        for(j in 1:N){
            M[i,j] <- (invlogit(a - d[i,j]) > .5)+0
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
##### Finding a good starting point #####


# hlsm.init <- function(){
#     list(alpha=c(1,1),
#          z=abind(b,b,along=3), # initializes the ovals as circles, also
#          b=b)
# }

library(expm) # has 'sqrtm'
procrustean <- function(z, z.0){
    Re(z.0 %*% t(z) %*% solve(sqrtm(z %*% t(z.0) %*% z.0 %*% t(z))) %*% z)
}

# x <- stretch(b, 2)
# plot.lsm(1, x)
# y <- stretch(b, 2, dim='y')
# plot.lsm(1, y, col=2, add=T)
# procrusted <- procrustean(x, y)
# plot.lsm(1, procrusted, col=2, add=T) # WORKS -- but is not the transform we need.


graph.error <- function(alpha, edges, z){
    # REMOVE THIS -- maybe it's fucking with optim?
#    if(is.null(alpha))
#        alpha <- 1
    g.hat <- positions.to.matrix(alpha, z)
    tbl <- table(as.vector(g.hat), as.vector(edges))
    # If the table collapses, then somethings's gone wrong. GODDAMMIT 100 AGAIN
    if(length(dim(tbl))!=2){
        return(100)
    } else {
        return(1-sum(diag(tbl))/sum(tbl))
    }
}

find.alpha <- function(g, z){
    return(optimize(partial(graph.error, edges=g, z=z),
                            lower=0, upper=5)$minimum)
}

# a <- find.alpha(two.ovals[,,1], b)
# 1-graph.accuracy(a, two.ovals[,,1], b) # THERE we go (but -- why not perfect? Oh, ovals.)
# plot.lsm(a, b) # what ... 

# # DO IT OUT THE FUCK MANUALLY
# out <- array(dim=100)
# v <- seq(0,50,length=100)
# for(i in 1:100){
#     out[i] <- graph.accuracy(v[i], two.ovals[,,1], b)
# }
# plot(out) # => 0.45. Okay ...
# optimize(partial(graph.accuracy, edges=two.ovals[,,1], z=b), lower=0, upper=5)$minimum
# plot.lsm(.45, b) # Okay, so THAT works. Good --




library(purrr) # partial()
# Do a separate lsm for each layer (making sure it)
# - has the right alpha => optim,
# - is the right scale for alpha=5
# THEN procrustean transform them all to match up with the center layer
find.init <- function(m){
    N <- dim(m)[1]
    K <- dim(m)[3]
    avg <- (apply(m, MARGIN=1:2, sum) > 0) + 0    # init failed with K=3; round=>0?
    lsm.avg <- ergmm(network(avg) ~ euclidean(d=2)^2) # Squared?
    z.avg <- lsm.avg$mcmc.mle$Z # .pmode v .mle? # 92 EDGES?
    # a.avg <- summary(lsm.avg)$pmean$coef.table[1]$Estimate # jesus christ, latentnet
     
    # NOW SCALE s.t. b has unit radius
    bounds <- apply(z.avg, 2, range)          # => x/y.min, x/y.max
    z.scale <- max(bounds[2,] - bounds[1,])/2 # => finds the max (AXIS ALIGNED) diameter
    z.avg <- z.avg * 1/z.scale                # (notice though, divide by HALF diam.)
    z.avg <- scale(z.avg, scale=F)            # and center
       
    z <- array(dim=c(N,K,2))    
    for(k in 1:K){
        lsm.k <- ergmm(network(m[,,k]) ~ euclidean(d=2)^2)
        z.k <- lsm.k$mcmc.mle$Z
        # a.k <- summary(lsm.k)$pmean$coef.table[1]$Estimate
        z.k <- anticrustean(z.k, z.avg) # Nobody expects anticrustes!
        z[,k,] <- z.k
    }
    
    # AND CHOOSE THE RIGHT ALPHA #
    # Now that we've RE-scaled
    
    alpha <- array(dim=(K)) # IS no alpha_b b/c no EDGES
    
    # alpha[1] <- find.alpha(avg, z.avg)
    for(k in 1:K){
        alpha[k] <- find.alpha(m[,,k], z[,k,])
    }
    
    init <- function(){
        return(list(alpha=alpha, z=z, b=z.avg, lsm=lsm.avg))
    }
    return(init)
}
# oval.init <- find.init(two.ovals)()
# plot.positions(oval.init$z, oval.init$b, oval.init$alpha) # 


# (Functions for ...)
##### RUNNING STAN #####

library(latentnet)


# Designed to go from parameters to model in one line
one.shot <- function(edges, init, sigma, iter, file, model_name){
    
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
    
    #hlsm.init <- find.init(edges)
    fit <- stan(file=file,
                model_name=model_name,
                iter=iter,
                data=hlsm.data,
                init=init, #hlsm.init,
                verbose=T,
                control=list(max_treedepth = 15))
    theta <- unpack.hlsm.fit(fit, N, K)#, which=which)
    return(list(fit=fit, theta=theta, init=init()))
}



# (Functions to ...)
##### EXAMINE THE ESTIMATED VALUES #####

unpack.hlsm.fit <- function(fit, N, K, which="max"){
    # s <- summary(fit)$summary[,"mean"]
    draws <- as.matrix(fit)
    if(which=="max")
        s <- draws[which.max(draws[,"lp__"]),]
    if(which=="min")
        s <- draws[which.min(draws[,"lp__"]),]
    if(which=="mean")
        s <- colMeans(draws)
    if(which=="last")
        s <- draws[nrow(draws),]
    
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


plot.positions <- function(z, b, alpha=NULL, xlim=NULL, ylim=NULL){
    # First, make sure the plot is big enough
    if(is.null(xlim))
        xlim <- range(c(z[,,1], b[,1]))
    if(is.null(ylim))
        ylim <- range(c(z[,,2], b[,2]))
    
    if(is.null(alpha)) alpha <- 5
    if(length(alpha)!=z[2])
        alpha <- rep(alpha, dim(z)[2]+1)
    
    plot.lsm(alpha[1], b, add=F, col=1, xlim=xlim, ylim=ylim)
    for(i in 1:dim(z)[2]){
        plot.lsm(alpha[i+1], z[,i,], add=T, col=i+1)
    }
}

plot.model <- function(m, alpha=NULL, which="max", xlim=NULL, ylim=NULL){
    
    # MOTHER FUCK
    x <- m$theta$z
    N <- dim(x)[1]
    K <- dim(x)[2]
    theta <- unpack.hlsm.fit(m$fit, N, K, which=which)
    
    z <- theta$z
    b <- theta$b
    init <- m$init$b
    
    if(is.null(alpha))
        alpha <- theta$alpha
    plot.positions(z, b, alpha=alpha, xlim=xlim, ylim=ylim)
    plot.lsm(alpha=1, z=init, col="blue", add=T)
}

####### FUNCTION TESTS #######
# plot.multigraph(zhat$z, alpha=.02)
# 
# plot.lsm(1, b)
# plot(b, xlim=c(-.5,.5), ylim=c(-.5,.5))
# plot.lsm(.02, zhat$b,      add=T)
# plot.lsm(.02, zhat$z[,1,], add=T, col="red") # .01 => disconnected
# plot.lsm(.02, zhat$z[,2,], add=T, col="blue")


# REMEMBER -- THIS MODEL IS STILL THE MEAN, NOT THE MAP (right?) 
model.accuracy <- function(edges, z, alpha=NULL){ # ooh, optim!
    K <- dim(edges)[3]
    if(is.null(alpha))
        alpha <- 1
    if(length(alpha)!=K)
        alpha <- rep(alpha, K)
    for(i in 1:dim(edges)[3]){
        g.hat <- positions.to.matrix(alpha[i], z[,i,])
        tbl <- table(as.vector(g.hat), as.vector(edges[,,i]))
        print(tbl)
    }
}

# model.accuracy(two.ovals, two.ovals.model, two.ovals.model$theta$alpha) # after ...
# model.accuracy(two.ovals, two.ovals.positions, 1) # okay so, perfect ... 
# model.accuracy(two.ovals, ovals.init$z, ovals.init$alpha)


# Draw little gray dotted lines between versions of the same node
plot.errors <- function(z, b){
    if(length(dim(z))==2){
        z <- array(z, dim=c(dim(z),1))
        z <- aperm(z, c(1,3,2))    
    }

    N <- dim(z)[1]
    K <- dim(z)[2]
    for(i in 1:N){
        for(k in 1:K){
            x <- c(z[i,k,1], b[i,1])
            y <- c(z[i,k,2], b[i,2])
            lines(x, y, col='gray', lty=2)
        }
    }
}
#plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)


#sigma <- function(x) 1/(1+exp(-x))
L22 <- function(x) sum(x^2) 

# gotta RE-run, becaue # gotta RE-run, because I changed the name
# of the parameters, z.hat => theta. Fuckin' ...
# model.accuracy(two.circles, two.circles.model)


##### SETTING PARAMETERS CREATING DATA AND RUNNING MODELS #####

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # => 4 chains at once. Issues?



##### INIT #####

N <- 20
circle <- seq(0,2*pi,length=N)
b <- cbind(cos(circle), sin(circle))
niter <- 1e3
sigma <- 10


library(abind)
library(igraph)
source('anticrustes.r') # necessary for find.init

# plot.lsm(1, b)
# plot.lsm(1, stretch(b, 2, dim="x"), add=T, col=2)
# plot.lsm(1, stretch(b, 2, dim="y"), add=T, col=4)
# plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)

# TWO ORTHOGONAL OVALS
two.ovals.positions <- aperm(abind(stretch(b, 2, dim="x"),
                                   stretch(b, 2, dim="y"), along=3),
                             c(1,3,2))
two.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="x"),
                                          stretch(b, 2, dim="y")))
plot.errors(two.ovals.positions, b) # Good!

#plot(graph_from_adjacency_matrix(two.ovals[,,1]))
#plot(graph_from_adjacency_matrix(two.ovals[,,2])) 

# SOME TESTS WITH -- REDUCING SIGMA. Issue, though, is -- how many iterations?
ovals.init <- find.init(two.ovals)() # done outside so you can plot them
plot.positions(ovals.init()$z, ovals.init()$b, ovals.init()$alpha)



##### FITTING #####
# two.ovals.model <- one.shot(two.ovals, ovals.init, .025, 4e4,
#                             "hlsm.stan", "two_ovals")
# plot.model(two.ovals.model)
# plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)





# AND NOW WITH THE TRUE PARAMETERS
# true.theta <- function()(list(a=c(1,1), b=b, z=two.ovals.positions))
# two.ovals.true.model <- one.shot(two.ovals, true.theta, .025, 2e4,
#                                  "hlsm.stan", "two_ovals_true")   
# plot.model(two.ovals.true.model, alpha=.5, which="max")
# plot.errors(two.ovals.true.model$theta$z, two.ovals.true.model$theta$b)
# 
# traceplot(two.ovals.model$fit) # I've never seen the z/e plots -- 
# model.accuracy(two.ovals, ovals.init()$z, ovals.init()$alpha) # Before: not BAD, 
# model.accuracy(two.ovals, two.ovals.model$theta$z, two.ovals.model$theta$alpha) # WORSE
# plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)

# PLOT THE INIT POSITIONS, TOO
# MAKE MOVES (coord-wise, x/y)
# VISUALIZE MOVES

























# 
# # plot.multigraph(two.ovals.model$theta$z)
# # z.hat <- unpack.hlsm.fit(three.ovals.model$fit, N, 3)
# 
# m <- two.ovals.model
# plot.positions(m$theta$z, m$theta$b, alpha=10)
# 
# # OKAY -- take a BUNCH, do them all ... except, fuckin ... fuck.
# # 'model' should be 'parameters,' no, 'positions,' plot 'fit' maybe? That's positions.
# plot.positions(two.ovals.positions, b, alpha=.1) # what the fuck ... 
# 
# 
# 
# # THREE ~ORTHOGONAL OVALS
# three.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="y"),
#                                             rotate(stretch(b, 2, dim="y"), pi/3),
#                                             rotate(stretch(b, 2, dim="y"),-pi/3)))
# three.ovals.model <- one.shot(three.ovals, sigma, "hlsm.stan", "three_ovals", niter)
# # 
# plot.model(three.ovals.model, alpha=c(.12, .1, .12), # ... blue?
#            xlim=c(-1,1), ylim=c(-1,1)) # oh VERY interesting -- what the hell?
# plot.multigraph(z.hat$z) # -- what
# plot.lsm(1, z.hat$z[,,1])
# 
# 
# # For the L1 test -- lambda?
# # TWO CONCENTRIC CIRCLES
# two.circles <- positions.to.multigraph(list(stretch(b, 2, dim=c("x", "y")),
#                                             b))
# two.circles.model <- one.shot(two.circles, sigma, "hlsm.stan", "two_circles", niter)
# 
# two.overlapping.circles   <- positions.to.multigraph(list(b, b))
# two.ovlp.circles.model   <- one.shot(two.overlapping.circles,  sigvec, "hlsm_lasso.stan", "two_ovlp_circles",    1e3)
# #three.ovlp.circles.model <- one.shot(three.ovlp.circles, sigvec, "hlsm_lasso.stan", "three_ovlp__circles", 2e4)
# 
# 
# three.overlapping.circles <- positions.to.multigraph(list(b, b, b))
# 
# 
# 
# 
# # What's ONE of these -- 
# K <- 2
# sigma <- 10
# sigmat <- array(0, dim=c(K,2,2))
# sigmat[,1,1] <- sigma
# sigmat[,2,2] <- sigma
# hlsm.data <- list(N=N,
#                   K=K,
#                   edges=edges,     # make sure these are ints  (not bool) else FLATTENED
#                   sigma_alpha=10,   # NO idea what these should be, using mvlsm's
#                   mu_b=c(0,0),
#                   sigma_b=sigma*diag(2), # [s   0; 0   s] (one matrix)
#                   sigma_z=sigmat)        # [s_k 0; 0 s_k] (one matrix per layer)
# 
# hlsm.data[["edges"]]   <- two.ovals
# hlsm.data[["sigma_z"]] <- sigmat
# hlsm.init <- find.init(two.ovals) # MAYBE THE PROBLEM ? TRY SIMPLER!
# plot.positions(hlsm.init()$z, b, 1) # This looks not too damn bad. What goes wrong?
# # How about -- plot over time? Do them one by one, watch the sampling progress?
# two.circles.fit <- stan(file="hlsm.stan",
#                         model_name="hlsm_sim_two_ovals",
#                         iter=2e4,
#                         data=hlsm.data, init=hlsm.init,
#                         verbose=T, control=list(max_treedepth = 15))
# # BEWARE: These are only defined BELOW. Will need to shuffle everything.
# zhat <- unpack.hlsm.fit(two.circles.fit, N, K, which="max")
# plot.positions(zhat$z, b, alpha=2)
# #plot.multigraph(zhat$z, zhat$b, alpha=5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### THE DUMP -- old tests, bad ideas, etc #####
# 
# plot.model(two.ovals.model, which="mean") # what the fuck ... what the FUCK though.
# plot.model(two.ovals.model, which="min") # the best (min?)
# plot.model(two.ovals.model, which="last") # the last
# 
# 
# # BARBELLS AND ROTATIONS -- maybe not very clear thinking
# # The sparse path connecting the dense ends
# w <- cbind(0.5*0:12-3, rep(0,13))
# # The two dense barbell ends
# x <- mvrnorm(10, mu=c(-3,0), Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
# y <- mvrnorm(10, mu=c(3,0),  Sigma=array(c(.1,0,0,.1), dim=c(2,2)))
# z <- rbind(w,x,y)
# plot(z)
# 
# 
# 
# # Finding a good starting point 0.1 #
# 
# # Better idea was to use latentnet itself -- MDS necessary, if you don't HAVE that.
# library(MASS) # thought that was automatic (loads isoMDS)
# # Sum the edges
# summed.edges <- apply(edges, MARGIN=1:2, sum) # THERE; also still symmetric
# d <- 1/summed.edges
# d[!is.finite(d)] <- 5 # Because. The biggest value in summed.edges = 5.
# mds <- isoMDS(as.matrix(d), k=2) # converged! That's good news -- 
# 
# hlsm.init <- function(){
#     # alpha, b, and z; using mds$points
#     z <- array(dim=c(N,K,2)) # defined above
#     for(k in 1:K)
#         z[,k,] <- lsm$mcmc.pmode$Z #mds$points
#     list(alpha=rnorm(K), z=z, b=lsm$mcmc.pmode$Z) #mds$points
# }
# 
# 
# 
# # ALL THIS WAS A GOOD FIRST TEST, BUT I DON'T NEED IT NOW --
# # Two ovals of data, at 90 degrees to one another. Hope to infer a circle of base positions
# edges <- array(dim=c(N,N,K))
# 
# # How do you stack 2d structures? ALONG=3 (below)
# # https://stackoverflow.com/questions/14857358/stacking-array-slices-rbind-for-2-dimensions
# library(abind)
# # Again, this isn't ... fuckin, exactly like the DGP, is it. Fuckin' ...
# positions <- abind(stretch(b, 3, dim="x"),
#                    stretch(b, 3, dim="y"), along = 3)
# for(i in 1:2){
#     edges[,,i] <- positions.to.matrix(positions[,,i], a=1)
# }
# 
# 
# 
# fit.20k <- stan(file="hlsm.stan",
#              model_name="hlsm_sim",
#              data=hlsm.data,
#              init=hlsm.init,
#              iter=2e4,
#              verbose=T,
#              control=list(max_treedepth = 15)) # from 10, as suggested
# # 2k was lousy, 5k better-looking (complete eyeball of the traceplot) 20k -- LOOKS ok
# # fit.2k <- fit
# # fit <- fit.20k
# # save(fit, file="hlsm_sim_fit_20k.Rdata")
# # zhat.table <- summary(fit)$summary[,"mean", drop=F]
# 
# # And finally, plot to check
# ghat.1 <- positions.to.matrix(.01, zhat$z[,1,])
# ghat.2 <- positions.to.matrix(.01, zhat$z[,2,])
# ghat.b <- positions.to.matrix(.01, zhat$b)
# 
# # plot(graph_from_adjacency_matrix(out)) # it's -- a big fat ball. Okay both are.
# # WAIT -- interesting -- if I just MAKE alpha=5, I get ~ the right answer
# # (if, remember, I also initialized it that way. Wait did I?)
# plot(graph_from_adjacency_matrix(positions.to.matrix(1, b)))
# 
# # ~Classification accuracies
# # What happened to 'edges'? It's not a table anymore
# table(as.vector(ghat.1), as.vector(edges[,,1])) # NOT TERRIBLE!
# table(as.vector(ghat.2), as.vector(edges[,,2])) # NOT TERRIBLE! ALSO!
# table(as.vector(ghat.b), as.vector(positions.to.matrix(b, a=1))) # Seriously not bad
# # Now try tuning with alpha, maybe? And why was it so high?
# # Also, try with initialization values that aren't cheating.
# 
# 
# out <- array(dim=20)
# vals <- seq(0, 5, length=20)
# for(i in 1:20){
#     cm <- table(as.vector(ghat.b), as.vector(positions.to.matrix(b, a=vals[i]))) 
#     out[i] <- sum(diag(cm))/sum(cm)
# }
# plot(out, xaxt="n")
# axis(1, at=1:20, labels=round(vals,2)) # Okay so yes, peaks around 1; why estimate=20?
# 
# 
