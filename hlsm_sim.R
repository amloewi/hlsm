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
rotate <- function(x, theta){
    f <- rotation.matrix(theta)
    return(f(x))
}
rotate.45 <- rotation.matrix(pi/4)
rotate.60 <- rotation.matrix(pi/3)
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
    #if(length(a)!=)
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

library(expm) # has 'sqrtm'
procrustean <- function(z, z.0){
    # REQUIRES CENTERING
    z   <- scale(z, scale=F)
    z.0 <- scale(z.0, scale=F)
    Re(z.0 %*% t(z) %*% solve(sqrtm(z %*% t(z.0) %*% z.0 %*% t(z))) %*% z)
}

# x <- stretch(b, 2)
# plot.lsm(1, x)
# y <- stretch(b, 2, dim='y')
# plot.lsm(1, y, col=2, add=T)
# procrusted <- procrustean(x, y)
# plot.lsm(1, procrusted, col=2, add=T) # THERE WE FUCKIN GO


graph.accuracy <- function(alpha, edges, z){
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
    return(optimize(partial(graph.accuracy, edges=g, z=z),
                            lower=0, upper=5)$minimum)
} # WORKS WHEN IT'S FIVE -- NOT ABOVE! WHAT THE WORTHLESS FUCK!
a <- find.alpha(two.ovals[,,1], b)
1-graph.accuracy(a, two.ovals[,,1], b) # THERE we go (but -- why not perfect? Oh, ovals.)
plot.lsm(a, b) # what ... 

# # DO IT OUT THE FUCK MANUALLY
# out <- array(dim=100)
# v <- seq(0,50,length=100)
# for(i in 1:100){
#     out[i] <- graph.accuracy(v[i], two.ovals[,,1], b)
# }
# plot(out) # => 0.45. Okay ...
# optimize(partial(graph.accuracy, edges=two.ovals[,,1], z=b), lower=0, upper=5)$minimum
# plot.lsm(.45, b) # Okay, so THAT works. Good --




library(purrr)
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
    a.avg <- summary(lsm.avg)$pmean$coef.table[1]$Estimate # jesus christ, latentnet
    # I should at least see if this matches with the -- optimal one.
    # Oh but, I'll need a new one anyway, if I -- re-scale.
        
    z <- array(dim=c(N,K,2))    
    for(k in 1:K){
        lsm.k <- ergmm(network(m[,,k]) ~ euclidean(d=2)^2)
        z.k <- lsm.k$mcmc.mle$Z # pmode v. mle? # 92 EDGES? 
        a.k <- summary(lsm.k)$pmean$coef.table[1]$Estimate # jesus christ, latentnet
        # Not quite done yet!
        z.k <- procrustean(z.k, z.avg) 
        z[,k,] <- z.k
    }
    
    # NOW SCALE
    # scale <- max(abs(range(z.avg)))
    bounds <- apply(z.avg, 2, range)
    scale <- max(bounds[2,] - bounds[1,])/2
    z.avg <- z.avg * 1/scale
    z     <- z     * 1/scale
    
    # AND CHOOSE THE RIGHT ALPHA #
    ##############################
    # Now that we've RE-scaled
    
    alpha <- array(dim=(K)) # IS no alpha_b b/c no EDGES
    
    # alpha[1] <- find.alpha(avg, z.avg)
    for(k in 1:K){
        alpha[k] <- find.alpha(m[,,k], z[,k,])
    }
    
    # FUCK -- if I don't have alpha anymore (or should I just fuckin' keep it)
    # THEN I -- need to tune the SCALE, to maximize representation. Hm. Fuck.

    init <- function(){
        # z <- array(dim=c(N,K,2))
        # for(k in 1:K)
        #     z[,k,] <- z.hat
        return(list(alpha=alpha, z=z, b=z.avg, lsm=lsm.avg))
    }
    return(init) # needs a FUNCTION, remember -- why, who knows.
}
oval.init <- find.init(two.ovals)()
plot.positions(oval.init$z, oval.init$b, oval.init$alpha)
# so -- optim returns 100 (the upper bound) but they all look real good at 1. WHY.
# One more reason to just say fuck alpha, no?

# (Functions for ...)
################
# RUNNING STAN #
################

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
################################
# EXAMINE THE ESTIMATED VALUES #
################################

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

# I NEED -- PLOT POSITIONS, plot.model calls that, and ... plot.multigraph
# uses igraph, and input is just the -- something. 

# This plots from a set of edges -- the graph itself. plot.model is for the estimate.
plot.multigraph <- function(Z, alpha=5){
    
    
    # plot.lsm(alpha, Z[,1,], add=F, col=1) # WATCH THE INDEXING HERE -- DID I CHANGE THIS?
    # for(i in 2:dim(Z)[3]){
    #     plot.lsm(1, Z[,2,], add=T, col=i)
    # }
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
model.accuracy(two.ovals, two.ovals.model, two.ovals.model$theta$alpha) # after ...
model.accuracy(two.ovals, two.ovals.positions, 1) # okay so, perfect ... 
model.accuracy(two.ovals, ovals.init$z, ovals.init$alpha)


# Draw little gray dotted lines between versions of the same node
plot.errors <- function(z, b){
    for(i in 1:N){
        for(k in 1:K){
            x <- c(z[i,k,1], b[i,1])
            y <- c(z[i,k,2], b[i,2])
            lines(x, y, col='gray', lty=2)
        }
    }
}
plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)

# gotta RE-run, becaue # gotta RE-run, because I changed the name
# of the parameters, z.hat => theta. Fuckin' ...
# model.accuracy(two.circles, two.circles.model)

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


# # # # # # # # # # # # # # # 
##############################
# RUNNING AND TESTING GRAPHS #
##############################

N <- 20
circle <- seq(0,2*pi,length=N)
b <- cbind(cos(circle), sin(circle))
niter <- 1e3
sigma <- 10


library(abind)
library(igraph)

# TWO ORTHOGONAL OVALS
plot.lsm(1, b)
plot.lsm(1, stretch(b, 2, dim="x"), add=T, col=2)
plot.lsm(1, stretch(b, 2, dim="y"), add=T, col=4)
plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)

two.ovals.positions <- aperm(abind(stretch(b, 2, dim="x"),
                                   stretch(b, 2, dim="y"), along=3),
                             c(1,3,2))
two.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="x"),
                                          stretch(b, 2, dim="y")))
plot.errors(two.ovals.positions, b) # Good!

plot(graph_from_adjacency_matrix(two.ovals[,,1]))
plot(graph_from_adjacency_matrix(two.ovals[,,2])) 
ovals.init <- find.init(two.ovals)

model.accuracy(two.ovals, ovals.init()$z, ovals.init()$alpha) # Before: not BAD, 
model.accuracy(two.ovals, two.ovals.model$theta$z, two.ovals.model$theta$alpha) # WORSE

# SOME TESTS WITH -- REDUCING SIGMA. Issue, though, is -- how many iterations?
ovals.init <- find.init(two.ovals) # done outside so you can plot them
plot.positions(ovals.init()$z, ovals.init()$b, ovals.init()$alpha)
two.ovals.model <- one.shot(two.ovals, ovals.init, .1, 2e4, # 4e3 not enough ... 2e4?
                            "hlsm.stan", "two_ovals")
traceplot(two.ovals.model$fit) # I've never seen the z/e plots -- 
# Also, gotta do better than the ... 
plot.model(two.ovals.model, which="max") # the best
# WHY ARE THEY ALL ON THAT FUCKIN DIAGONAL? Am I even DOING this right?
plot.errors(two.ovals.model$theta$z, two.ovals.model$theta$b)

# Hm ... LOOKS good, but the errors are shitty. Hm ... procrustes might be the wrong
# transpose. I don't want to match the shape, I want to match the node positions --
# is it doing that? How do I minimize the distances between THEM?

# PLOT THE INIT POSITIONS, TOO
# MAKE MOVES (coord-wise, x/y)
# VISUAL MOVES

my.transform <- function(){
    # Rotate, until you min distances between original points
    # translate, for same
    # Is this the same as procrustes? Should look it up --
    # ALSO -- coordinate descent/convexity for LSMs in general
}


# uniroot! What a great name!

sigma <- function(x) 1/(1+exp(-x))
L22 <- function(x) sum(x^2) 

# For speed, I should check -- this against other things
build.b.fxn <- function(y, a, b, e, i, axis){
    f <- function(x){
        total <- 0
        for(j in (1:N)[-i]){ # Does this even matter? As long as the intercept's consistent
            for(k in 1:K){
                bi <- b[i,]
                bi[axis] <- x
                
                d <- bi + e[i,k,] - b[j,] - e[j,k,]
                r <- sigma(a - L22(d)) - y[i,j,k] # This should be a scalar
                
                total <- total + r*sqrt(L22(d)) # right, THIS. Fuck. Kay ... d is the problem
                
            }
        }
        return(-total)
    }
    return(f)
}

build.e.fxn <- function(y, a, b, e, i, k, axis, l){
    f <- function(x){
        total = 0
        for(j in (1:N)[-i]){ # Does this even matter? As long as the intercept's consistent
            eik <- e[i,k,]
            eik[axis] <- x
            d <- b[i,] + eik - b[j,] - e[j,k,]
            r <- sigma(a - L22(d)) - y[i,j,k]
            total <- total + r*sqrt(L22(d)) # + -l/l #? HOW DO I SOLVE THESE? manually? line search?
        }
        return(-total)
    }
    return(f)
}

lkhd <- function(y, a, b, e, l){
    total <- 0
    N <- dim(y)[1]
    K <- dim(y)[3]
    for(k in 1:K){
        for(i in 1:N){
            for(j in 1:N){
                s <- sigma(a[k] - L22(b[i,] + e[i,k,] - b[j,] - e[j,k,]))
                y. <- y[i,j,k]    
                total <- total + s^y.*(1-s)^(1-y.)
            }
        }
    }
    total <- total - l*sum(abs(e)) # right -- subtract, for max (REMEMBER IT'S MAX)
}

coord.opt <- function(y, a, b, e, l, maxit=1){

    iter <- 1
    lkhd <- array(dim=maxit)
    while(iter < maxit){ # Also a convergence criterion -- what? 
        # For each layer, first update the layer, then the base
        for(k in 1:K){ 
            # update the e_ik
            for(i in 1:N){
                for(axis in 1:2){
                    f.bi <- build.b.fxn(y, a[k], b, e, i, axis)
                    #bi.star <- uniroot(f.bi, c(0,.), extendInt = "yes")$root
                    bi.star <- optimize(f.bi, c(-20,20))$minimum
                    b[i,axis] <- bi.star
                }
            }
            # update the b_i
            for(i in 1:N){
                for(axis in 1:2){
                    f.eik <- build.e.fxn(y, a[k], b, e, i, k, l, axis)
                    # eik.star <- uniroot(f.eik, c(9,10), extendInt = "downX")$root
                    # eik.star <- optimize(f.eik, c(-20,20))$minimum
                    
                    obj.min <- Inf
                    eik.star <- NULL
                    for(v in seq(-20,20,length=1000)){
                        obj.val <- f.eik(v)
                        if(obj.val < obj.min){
                            obj.min <- obj.val
                            eik.star <- v
                        }
                        if(obj.val > -l & obj.val < l){ # if v \in [-lambda, lambda]
                            eik.star <- 0
                            break
                        }
                    }
                    e[i,k,axis] <- eik.star
                    
                }
            }
        }
        lkhd[iter] <- lkhd(y, a, b, e, l)
        iter <- iter + 1
    }
    theta <- list(a=a, b=b, e=e, l=l, lkhd=lkhd)
    return(theta)
}


# THERE IS no 'e' variable right now -- e = z - b, but the indices are a little tricky
a  <- ovals.init()$alpha
b. <- ovals.init()$b # CAREFUL -- don't overwrite the REAL b
z  <- ovals.init()$z
e  <- array(dim=c(N,K,2))
for(k in 1:K){
    e[,k,] <- z[,k,] - b.
}
lambda <- 1
f.bi  <- build.b.fxn(two.ovals, 1, b, e, 1, 1)
f.eik <- build.e.fxn(two.ovals, 1, b, e, 1, 1, 1, 1) # fuck, can be maximixed (~11)

out <- optimize(f.bi, c(0,20)) # except THAT was ok -- hm, got the edge though -- 

plot(sapply(-10:10, f.bi)) # hm -- negative forever? Oh also, MAX OR MIN?
# ... isn't what what I'm trying do -- do what with, now?
ur <- uniroot(f.bi, c(-3,0), extendInt = "yes") # WHAT THE FUCK, UNIROOT
ur$root


abline(h=0)
out$minimum # what ... .47 ? oh it does it by INDEX .. wait now ... oh btwn 1/2. Yes!

f.eik(-1) # ah, makes a VECTOR -- hm. f.bi(0) # ah, makes a VECTOR -- hm. 

source("/path/to/file/my_fn_lib1.r")


# I should really check -- all sorts of things. lkhd, accuracy, effectiveness
# of the individual optimizations ... BUT, it looks like it's RUNNING. At least.

# with lambda=1, everything got zeroed out. Dropping to .1.
out <- coord.opt(two.ovals, a, b., e, .1, maxit=20) # FUCK. What is it and why? 'a'.
# mother fuck -- it appears to have ... is the penalty on backwards? Back up, asshole.


# Okay well now, it's all OVER the fucking place. K -- fuck. Need to back up,
# commit, and try again. Fuck.
plot(out$lkhd) # IT IS INCREASING. LET'S SEE HOW LONG. (Does it plateau?)

# Hm -- peaked soon, then -- dropped precipitously. Hm. K. What was happening?
# WHAT'S HAPPENING WITH ACCURACY?
plot.positions(ovals.init()$z, ovals.init()$b, ovals.init()$alpha)
z.out  <- array(dim=c(N,K,2))
for(k in 1:K){
    z.out[,k,] <- out$b + out$e[,k,]
}
plot.positions(z.out, out$b, .4) # why is alpha still a fucking problem?
plot.errors(z.out, out$b) 










# Okay so that doesn't work -- since it's 1d, do some plots of the f.bi, f.eik etc.,
# see if they HAVE roots (among other things)









# plot.multigraph(two.ovals.model$theta$z)
# z.hat <- unpack.hlsm.fit(three.ovals.model$fit, N, 3)

m <- two.ovals.model
plot.positions(m$theta$z, m$theta$b, alpha=10)

# OKAY -- take a BUNCH, do them all ... except, fuckin ... fuck.
# 'model' should be 'parameters,' no, 'positions,' plot 'fit' maybe? That's positions.
plot.positions(two.ovals.positions, b, alpha=.1) # what the fuck ... 








# THREE ~ORTHOGONAL OVALS
three.ovals <- positions.to.multigraph(list(stretch(b, 2, dim="y"),
                                            rotate(stretch(b, 2, dim="y"), pi/3),
                                            rotate(stretch(b, 2, dim="y"),-pi/3)))
three.ovals.model <- one.shot(three.ovals, sigma, "hlsm.stan", "three_ovals", niter)
# 
plot.model(three.ovals.model, alpha=c(.12, .1, .12), # ... blue?
           xlim=c(-1,1), ylim=c(-1,1)) # oh VERY interesting -- what the hell?
plot.multigraph(z.hat$z) # -- what
plot.lsm(1, z.hat$z[,,1])






# For the L1 test -- lambda?
# TWO CONCENTRIC CIRCLES
two.circles <- positions.to.multigraph(list(stretch(b, 2, dim=c("x", "y")),
                                            b))
two.circles.model <- one.shot(two.circles, sigma, "hlsm.stan", "two_circles", niter)

two.overlapping.circles   <- positions.to.multigraph(list(b, b))
two.ovlp.circles.model   <- one.shot(two.overlapping.circles,  sigvec, "hlsm_lasso.stan", "two_ovlp_circles",    1e3)
#three.ovlp.circles.model <- one.shot(three.ovlp.circles, sigvec, "hlsm_lasso.stan", "three_ovlp__circles", 2e4)


three.overlapping.circles <- positions.to.multigraph(list(b, b, b))








# What's ONE of these -- 
K <- 2
sigma <- 10
sigmat <- array(0, dim=c(K,2,2))
sigmat[,1,1] <- sigma
sigmat[,2,2] <- sigma
hlsm.data <- list(N=N,
                  K=K,
                  edges=edges,     # make sure these are ints  (not bool) else FLATTENED
                  sigma_alpha=10,   # NO idea what these should be, using mvlsm's
                  mu_b=c(0,0),
                  sigma_b=sigma*diag(2), # [s   0; 0   s] (one matrix)
                  sigma_z=sigmat)        # [s_k 0; 0 s_k] (one matrix per layer)

hlsm.data[["edges"]]   <- two.ovals
hlsm.data[["sigma_z"]] <- sigmat
hlsm.init <- find.init(two.ovals) # MAYBE THE PROBLEM ? TRY SIMPLER!
plot.positions(hlsm.init()$z, b, 1) # This looks not too damn bad. What goes wrong?
# How about -- plot over time? Do them one by one, watch the sampling progress?
two.circles.fit <- stan(file="hlsm.stan",
                        model_name="hlsm_sim_two_ovals",
                        iter=2e4,
                        data=hlsm.data, init=hlsm.init,
                        verbose=T, control=list(max_treedepth = 15))
# BEWARE: These are only defined BELOW. Will need to shuffle everything.
zhat <- unpack.hlsm.fit(two.circles.fit, N, K, which="max")
plot.positions(zhat$z, b, alpha=2)
#plot.multigraph(zhat$z, zhat$b, alpha=5)









#########################################
#########################################
# THE DUMP -- old tests, bad ideas, etc #
#########################################
#########################################

plot.model(two.ovals.model, which="mean") # what the fuck ... what the FUCK though.
plot.model(two.ovals.model, which="min") # the best (min?)
plot.model(two.ovals.model, which="last") # the last


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
# WAIT -- interesting -- if I just MAKE alpha=5, I get ~ the right answer
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


