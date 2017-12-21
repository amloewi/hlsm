####################################
# REQUIRES FUNCTIONS IN hlsm_sim.R #
####################################


invlogit <- function(x) 1/(1+exp(-x))
L22 <- function(x) sum(x^2) 

# For speed, I should check -- this against other things
build.b.fxn <- function(y, a, b, e, i, axis){
    f <- function(x){
        total <- 0
        
        bi <- b[i,]
        bi[axis] <- x
        for(j in (1:N)[-i]){ # Does this even matter? As long as the intercept's consistent
            for(k in 1:K){
                d <- bi + e[i,k,] - b[j,] - e[j,k,]
                r <- invlogit(a - L22(d)) - y[i,j,k] # This should be a scalar
                
                total <- total + r*sqrt(L22(d))
            }
        }
        return(total)
    }
    return(f)
}

build.e.fxn <- function(y, a, b, e, i, k, axis, l){
    f <- function(x){
        total = 0
        
        eik <- e[i,k,]
        eik[axis] <- x
        for(j in (1:N)[-i]){
            d <- b[i,] + eik - b[j,] - e[j,k,]
            r <- invlogit(a - L22(d)) - y[i,j,k]
            
            total <- total + r*sqrt(L22(d))
        }
        return(total) # add the interval +/-lambda when trying to solve, outside
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
                s <- invlogit(a[k] - L22(b[i,] + e[i,k,] - b[j,] - e[j,k,]))
                y. <- y[i,j,k]    
                total <- total + s^y.*(1-s)^(1-y.)
            }
        }
    }
    total <- total - l*sum(abs(e)) # right -- subtract, for max (REMEMBER IT'S MAX)
    return(total)
}



f.argmax.b <- function(y, a, b, e, l, i, axis){
    f <- function(x){
        # make it just a function of -- ? then do one for eik (SO TIRED)
        b[i,axis] <- x
        return(lkhd(y, a, b, e, l))
    }
    return(f)
}

argmax.b <- function(y, a, b, e, l, i, axis){
    bi.old <- b[i,axis]
    f <- f.argmin.b(y, a, b, e, l, i, axis)
    v <- seq(-5, 5, length=1000) 
    f.v <- sapply(v, f)
    bi.star <- v[which.max(f.v)]
    if(is.null(bi.star)) bi.star <- bi.old
    return(bi.star)
}


f.argmax.e <- function(y, a, b, e, l, i, k, axis){
    f <- function(x){
        # make it just a function of -- ? then do one for eik (SO TIRED)
        e[i,k,axis] <- x
        return(lkhd(y, a, b, e, l))
    }
    return(f)
}

argmax.e <- function(y, a, b, e, l, i, k, axis){
    eik.old <- e[i,k,axis]
    f <- f.argmin.e(y, a, b, e, l, i, k, axis)
    v <- seq(-5, 5, length=1000) 
    f.v <- sapply(v, f)
    eik.star <- v[which.max(f.v)]
    if(is.null(eik.star)) eik.star <- eik.old
    return(eik.star)
}



cow <- function(y, a, b, e, l, maxit=2){
    lkhds <- array(dim=maxit)
    for(iter in 1:maxit){
        for(k in 1:K){                                      # For each layer, 
            print(paste("Layer: ", k))
            for(i in 1:N){
                print(paste("Node: ", i))
                for(axis in 1:2){
                    b[i,axis] <- argmax.b(y, a, b, e, l, i, axis)
                }
            }
            for(i in 1:N){
                for(axis in 1:2){
                    e[i,k,axis] <- argmax.e(y, a, b, e, l, i, k, axis)
                }
            }
        }
        lkhds[iter] <- lkhd(y, a, b, e, l)
    }
    return(list(b=b, e=e, lkhds=lkhds)) # And I still haven't opdated 'a'
}

# cowt <- cor(two.ovals, a, b., e, l, 2)





coord.opt <- function(y, a, b, e, l, maxit=1, trace=F){
    
    b.mle <- array(dim=dim(b))
    e.mle <- array(dim=dim(e))
    mle <- -Inf
    
    if(trace){
        pdf("coord_opt_trace.pdf")
        par(mfrow=c(3,2))
    }
    
    iter <- 1
    lkhd <- array(dim=maxit)
    while(iter < maxit){ # Also a convergence criterion -- what though?
        print(paste("Iter #:", iter))
        # For each layer, first update the base, then the layer 
        for(k in 1:K){ 
            # update the b_i
            for(i in 1:N){
                for(axis in 1:2){
                    f.bi <- build.b.fxn(y, a[k], b, e, i, axis)
                    
                    # bi.star <- broot.finding(f.bi, -10, 10)
                    # Make a heat-map overlay, of the values?
                    #print(paste("bi.star: ", bi.star))
                    
                    v <- seq(-5, 5, length=1000)
                    f.bi.v <- sapply(v, f.bi)
                    bi.star <- v[which.min(abs(f.bi.v))] # closest to 0 as possible

                    if(is.null(bi.star)) bi.star <- b[i,axis] # don't change it? else?
                    
                    if(trace){
                        snapshot(e, b, a, i, axis, bi.star)
                        
                        plot(f.bi.v, xaxt='n', main="b") # the values being chosen by brootfind
                        axis(1, at=100*(1:5)+1, labels=round(v[100*(1:5)+1], 1))
                        abline(v=which.min(abs(f.bi.v))) # i.e. the INDEX of the place
                    }
                    
                    b[i,axis] <- bi.star
                }
            }
            # update the e_ik
            for(i in 1:N){
                for(axis in 1:2){
                    f.eik <- build.e.fxn(y, a[k], b, e, i, k, l, axis)
                    # eik.star <- uniroot(f.eik, c(9,10), extendInt = "downX")$root
                    # eik.star <- optimize(f.eik, c(-20,20))$minimum
                    
                    # eik.star <- broot.finding(f.eik, -5, 5, l)
                    # print(paste("eik.star: ", eik.star))

                    v <- seq(-5, 5, length=1000)
                    f.eik.v <- sapply(v, f.eik)
                    eik.star <- v[which.min(abs(f.eik.v))] # closest to 0 as possible
                    
                    if(is.null(eik.star)) eik.star <- e[i,k,axis] # don't change it?
                    
                    if(trace){
                        snapshot(e, b, a, i, axis, eik.star, k)
                        
                        plot(f.eik.v, xaxt='n', main="e") # the values being chosen by brootfind
                        axis(1, at=100*(1:5)+1, labels=round(v[100*(1:5)+1], 1))
                        abline(v=which.min(abs(f.eik.v)))   # i.e. the INDEX of the place                    
                        abline(v=which.min(abs(f.eik.v-l)), col='blue') # that plus the boundary
                        abline(v=which.min(abs(f.eik.v+l)), col='red') # that minues the boundary
                    }
                    e[i,k,axis] <- eik.star
                    
                }
            }
        }
        
        # DO A MORE FINE-GRAINED LKHD GRAPH -- step by step, not just iteration. Right.
        
        lkhd[iter] <- lkhd(y, a, b, e, l)
        if(lkhd[iter] > mle){
            mle <- lkhd[iter]
            b.mle <- b
            e.mle <- e
        }
        print(paste('Lkhd: ', lkhd[iter])) # oh also, at each step? Maybe eq'ns wrong?
        
        iter <- iter + 1
    }
    if(trace)
        dev.off()
    
    theta <- list(a=a, b=b.mle, e=e.mle, l=l, lkhd=lkhd)
    return(theta)
}

broot.finding <- function(f, lower, upper, tol=.15, length=1000){
    v <- seq(lower, upper, length=length)
    # f.v <- sapply(v, f) # Unnecessary b/c short-circuiting -- find A solution.
    # If there are several, don't know which to choose anyway.
    for(i in 1:length){
        f.v <- f(v[i])
        if(-tol < f.v & f.v < tol){
            return(v[i])
        }
    }
    # Else what?
}

e.to.z <- function(e, b){
    
    N <- dim(e)[1]
    K <- dim(e)[2]
    
    z.out  <- array(dim=c(N,K,2))
    for(k in 1:K){
        z.out[,k,] <- b + e[,k,]
    }    
    return(z.out)
}

z.to.e <- function(z, b){
    
    N <- dim(z)[1]
    K <- dim(z)[2]
    
    e.out  <- array(dim=c(N,K,2))
    for(k in 1:K){
        e.out[,k,] <- z[,k,] - b
    }    
    return(e.out)
}

# Run WITH theta.star, but before SETTING it; allows for before/after
snapshot <- function(e, b, a, i, axis, star, k=NULL){
    z <- e.to.z(e, b)
    plot.positions(z, b, a) # add=F
    if(is.null(k)){
        bi.old <- b[i,]
        points(b[i,1], b[i,2], pch=19, col=4)
        
        b[i,axis] <- star
        plot.positions(z, b, a) # NEW graph; two per line/move, before+after
        points(bi.old[1], bi.old[2], pch=1, col=4)
        points(b[i,1], b[i,2], pch=19, col=4)
        lines(c(bi.old[1], b[i,1]),
              c(bi.old[2], b[i,2]), lty=2, col=4)
    } else {
        eik.old <- e[i,k,]
        points(e[i,k,1], e[i,k,2], pch=19, col=4)
        
        e[i,k,axis] <- star
        plot.positions(z, b, a) # NEW graph; two per line/move, before+after
        points(eik.old[1], eik.old[2], pch=1, col=4)
        points(e[i,k,1], e[i,k,2], pch=19, col=4)
        lines(c(eik.old[1], e[i,k,1]),
              c(eik.old[2], e[i,k,2]), lty=2, col=4)
    }
}


# bi.star <- optimize(f.bi, c(-20,20))$minimum






# THERE IS no 'e' variable right now -- e = z - b, but the indices are a little tricky
N <- 20
a  <- ovals.init()$alpha
b. <- ovals.init()$b # CAREFUL -- don't overwrite the REAL b
z  <- ovals.init()$z
K <- 2
e  <- array(dim=c(N,K,2))
for(k in 1:K){
    e[,k,] <- z[,k,] - b.
}
lambda <- 1

#############################################
cowt <- cow(two.ovals, a, b., e, lambda, 10) # SHE RUNS. Glacial. Look into it tomorrow.
# 12:14:44 ... ~ 10. seconds per 'Node' loop. ~=800s/it. Run a bunch, go to sleep!
#############################################
plot.positions(e.to.z(cowt$e, cowt$b), cowt$b) # WHAT the FUCK




setwd("~/Documents/Research/gsa/study/data")
load('multigraph.Rdata')
multigraph.init <- find.init(Y)
gsa.1 <- run.optimization(multigraph.init()$b,
                        aperm(e.to.z(multigraph.init()$z,
                                     multigraph.init()$b),
                              c(1,3,2)),
                        Y,
                        lambda=1,
                        t=1e-3, 
                        group.level=FALSE, n.steps=1500)
plot.positions(e.to.z(aperm(gsa.1$epsilon, c(1,3,2)), gsa.1$b), gsa.1$b, alpha=.1,
               xlim=c(-1,1), ylim=c(-1,1), one.plot=F)

# I WANT TO DO THIS FOR -- 1) MANY LAMBDA VALUES, 2) WITH LAYERS SEPARATED 3) w/ alpha.hats
# CONVERGENCE:
# lambda=1: 1500 iterations
gsa.01 <- run.optimization(multigraph.init()$b,
                        aperm(e.to.z(multigraph.init()$z,
                                     multigraph.init()$b),
                              c(1,3,2)),
                        Y,
                        lambda=.1,
                        t=1e-3, 
                        group.level=FALSE, n.steps=3000)
plot.positions(e.to.z(aperm(gsa.01$epsilon, c(1,3,2)), gsa.01$b), gsa.01$b, alpha=.1,
               one.plot=F, xlim=c(-1,1), ylim=c(-1,1))

gsa.10 <- run.optimization(multigraph.init()$b,
                           aperm(e.to.z(multigraph.init()$z,
                                        multigraph.init()$b),
                                 c(1,3,2)),
                           Y,
                           lambda=10,
                           t=1e-3, 
                           group.level=FALSE, n.steps=3000)
plot.positions(e.to.z(aperm(gsa.10$epsilon, c(1,3,2)), gsa.10$b), gsa.10$b, alpha=.1,
               one.plot=F)

gsa.100 <- run.optimization(multigraph.init()$b,
                           aperm(e.to.z(multigraph.init()$z,
                                        multigraph.init()$b),
                                 c(1,3,2)),
                           Y,
                           lambda=100,
                           t=1e-3, 
                           group.level=FALSE, n.steps=500)
plot.positions(e.to.z(aperm(gsa.100$epsilon, c(1,3,2)), gsa.100$b), gsa.100$b, alpha=.1)

                             



out <- coord.opt(two.ovals, a, b., e, 1, maxit=10, trace=F)


plot(out$lkhd) # up or down? up, right? Up, good -- but RESULTS. 
z.out  <- array(dim=c(N,K,2))
for(k in 1:K){
    z.out[,k,] <- out$b + out$e[,k,]
}
plot.positions(z.out, out$b, 1) # COMPLETELY fucked. (Also -- )
plot.errors(z.out, out$b)


# MAKE A -- GRAPH OF THE CHANGES. AS IT MOVES. pdf inside the loops, snapshots.



f.bi  <- build.b.fxn(two.ovals, 1, b, e, 1, 1) # => MAX/MIN ... what SHOULD it be?
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
plot.positions(z.out, out$b, .4) # COMPLETELY fucked
plot.errors(z.out, out$b) 



# THIS SHOULD BE ROOT FINDING BUT ITS NOT WORKING

# FIND INTERVAL? With a line search? (Just -- use the line search?)
# ALSO, these will be a lot easier to search with scaled versions.
# SO -- work on scaling. Fuck. 
# 
# tryCatch({
#     obj.vals <- array(dim=1000)
#     vals <- seq(-20,20,length=1000)
#     for(v.i in 1:1000){
#         obj.vals[v.i] <- f.bi(vals[v.i])
#     }
#     lower <- vals[which.min(obj.vals)]
#     upper <- vals[which.max(obj.vals)]
#     print(lower)
#     print(upper)
#     bi.star <- uniroot(f.bi, c(lower,upper), extendInt = "yes")$root#bi.star <- uniroot(f.bi, c(0,.), extendInt = "yes")$root#bi.star <- uniroot(f.bi, c(0,.), extendInt = "yes")$root
# })

# eik.star <- 0
# for(v in seq(-20,10,length=1000)){
#     val <- f.eik(v)
#     if((-lambda < val) & (val < lambda)){
#         eik.star <- v
#         break
#     }
# }

# obj.min <- Inf
# eik.star <- NULL
# for(v in seq(-20,20,length=1000)){
#     obj.val <- f.eik(v)
#     if(obj.val < obj.min){
#         obj.min <- obj.val
#         eik.star <- v
#     }
#     if(obj.val > -l & obj.val < l){ # if v \in [-lambda, lambda]
#         eik.star <- 0
#         break
#     }
# }
