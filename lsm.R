setwd("~/Documents/Classes/10-725/hlsm")
# GENERATING LSM DATA #
# The sparse path connecting the dense ends
w <- cbind(0.5*0:8-2, rep(0,9))
# The two dense barbell ends
x <- mvrnorm(10, mu=c(-2,0), Sigma=array(c(.2,0,0,.2),dim=c(2,2)))
y <- mvrnorm(10, mu=c(2,0), Sigma=array(c(.2,0,0,.2),dim=c(2,2)))
z <- rbind(w,x,y)

d <- as.matrix(dist(z))
a <- 1
N <- nrow(z)
X <- array(dim=c(N,N))

logit <- function(x) 1/(1+exp(-x))

for(i in 1:N){
    for(j in 1:N){
        # This is what a LSM assumes
        X[i,j] <- (logit(a - d[i,j]) > .5)+0
    }
}
mean(X)

lsm.data <- list(N=N,
                 edges=X,
                 sigma_alpha=10,
                 mu_z=c(0,0),
                 sigma_z=array(c(20,0,0,20), dim=c(2,2)),
                 sigma_fixed=array(c(.5,0,0,.5), dim=c(2,2)))

# This initialization provides the TRUE values of the parameters
lsm.init <- function(){
    list(alpha=a, z=z)
}

library(rstan)
lsm.fit <- stan("lsm.stan",
                model_name="lsm",
                data=lsm.data,
                init=lsm.init,
                verbose=T)

s <- summary(lsm.fit)$summary[,1]
a.hat <- s[1]
s <- s[2:(length(s)-1)] # cut off 'alpha' AND 'lp__'
x.hat <- s[c(T,F)]
y.hat <- s[c(F,T)]
z.hat <- cbind(x.hat, y.hat)
plot(y.hat ~ x.hat, z.hat, # z[,2] ~ z[,1], #
     col=c(rep('blue', 9), rep('red', 10), rep('green', 10)))
plot(z[,2] ~ z[,1], #
     col=c(rep('blue', 9), rep('red', 10), rep('green', 10)))


# Here is an independent verification of the data -- this model fits fine,
# and a plot of the result looks exactly as I think it should (a barbell).
#library(latentnet)
#lsm.ergmm <- ergmm(X ~ euclidean(d=2))
#plot(lsm.ergmm)
