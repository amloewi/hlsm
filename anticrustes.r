dilate <- function(x, a){
    return((x %*% diag(2))*a)
}
plot.lsm(1, dilate(b, 2))
plot.lsm(1, dilate(b, .5), add=T)
plot.errors(dilate(b,2), dilate(b,.5)) # LOOKS accurate enough. 


# 'a' is a length-2 vector, x is Nx2. 
translate <- function(x, a){
    return(t(t(x)+a))
}

translate.x <- function(x, a){
    x[,1] <- x[,1] + a
    return(x)
}

translate.y <- function(x, a){
    x[,2] <- x[,2] + a
    return(x)
}

# These two, pasted from _sim
rotation.matrix <- function(theta){
    R <- array(c(cos(theta),  sin(theta),
                 -sin(theta), cos(theta)), dim=c(2,2))
    return(function(x){x%*%R})
}
rotate <- function(x, theta){
    f <- rotation.matrix(theta)
    return(f(x))
}

posn.diff <- function(z, b){
    total <- 0
    for(i in 1:N){
        x <- z[i,] #c(z[i,1], b[i,1])
        y <- b[i,] #c(z[i,2], b[i,2])
        #lines(x, y, col='gray', lty=2)
        total <- total + L22(x - y)
    }
    return(total)
}

transformer <- function(z, b, f){
    g <- function(v){
        posn.diff(f(z, v), b)
    }
    return(g)
}

b.dilate <- b*2
dilate.star <- optimize(transformer(b.dilate, b, dilate), c(.1,5))$minimum # .1; bad.
plot.lsm(1, b)
plot.lsm(1, b.dilate, add=T)
plot.lsm(1, dilate(b.dilate, dilate.star), col=2, add=T) # THERE we go


b.rotate <- rotate.45(stretch(b,2))
rotate.star <- optimize(transformer(b.rotate, stretch(b,2), rotate), c(0,2*pi))$minimum
# 5.5; ?
plot.lsm(1, stretch(b,2))
plot.lsm(1, b.rotate, col=2, add=T)
plot.lsm(1, rotate(b.rotate, rotate.star), add=T, col=4) # YES!


b.x <- stretch(b, 2, "x")
b.y <- stretch(b, 2, "y")
plot.lsm(1, b.x)
plot.lsm(1, b.y, add=T, col=2)
plot.errors(b.x, b.y) # Right. Now, we see if this is minimizing.
plot(sapply(seq(0, 2*pi, length=360), transformer(b.y, b.x, rotate))) # looks good.

b.x.rotate <- rotate.45(b.x) # so we fuck it up ... see if it returns. Right.

rotate.x.star <- optimize(transformer(b.x.rotate, b.y, rotate), c(0,2*pi))$minimum
plot.lsm(1, b.y)
plot.lsm(1, rotate(b.x.rotate, rotate.x.star), add=T, col=4) # YES!


b.translate <- b+3
translate.star <- optimize(transformer(b.translate, b, translate.x), c(-5,5))$minimum
plot.lsm(1, b)
plot.lsm(1, b.translate, col=2, add=T)
plot.lsm(1, translate.x(b.translate, translate.star), add=T, col=4)

# OKAY THEN -- wrap them all together. (And maybe test against ... procrustes. 
# Could I have just maybe just read the paper? No, I think I -- found it's a pretty
# different transform. Good.)

anticrustean <- function(z, z.0){
    
    # Rotation needs to be first, AND it has to be APPLIED, before translations
    r.star   <- optimize(transformer(z, z.0, rotate),      c(  0,2*pi))$minimum
    z <- rotate(z, r.star)
    
#    for(i in 1:1){ # This goes ALL to hell when i > 1, but -- it doesn't look DONE.
    t.x.star <- optimize(transformer(z, z.0, translate.x), c(-20,  20))$minimum
    t.y.star <- optimize(transformer(z, z.0, translate.y), c(-20,  20))$minimum
    z <- translate.x(z, t.x.star)
    z <- translate.y(z, t.y.star)

    # This one has been weird -- maybe now that the order is right? Looks good!
    # d.star   <- optimize(transformer(z, z.0, dilate),      c(0.01,  3))$minimum
    d.star <- brute.optimize(transformer(z, z.0, dilate), 0.1,  10)
    print(d.star)
    z <- dilate(z, d.star)

    return(z)
}

# Fuckin -- fine, test the function overall. Like ...
plot.lsm(5, oval.init$z[,1,], xlim=c(-5,10), ylim=c(-10,0))
plot.lsm(5, oval.init$z[,2,], xlim=c(-5,10), ylim=c(-10,0),   add=T, col=2)
plot.lsm(5, anticrustean(oval.init$z[,1,], oval.init$z[,2,]), add=T, col=3)
plot.errors(anticrustean(oval.init$z[,1,], oval.init$z[,2,]), oval.init$z[,2,])

plot(sapply(seq(.01, 3, length=30),
            transformer(ovals.init()$z[,1,], ovals.init()$b, dilate)))
seq(.01, 3, length=30)[4] # .42 -- so ... (This is not working.)
# So do a -- line search? My own, fucking, fuck this function, version.

brute.optimize <- function(f, lower, upper, length=1000){
    vals <- seq(lower, upper, length=length)
    vals[which.min(sapply(vals, f))]
}
brute.optimize(transformer(ovals.init()$z[,1,], ovals.init()$b, dilate), .01, 10)
# => .25. What's the answer though? Well just test it.

# I THINK the transform is good -- problem is the search RANGES. Am I pinning alpha?
# And how does this play with ... alpha searches? Oh right, if we CAN adjust it ...
# 1) scale b to unit radius (and center)
# 2) anticrustes all the layers, to transformed b
# 3) find appropriate alphas for b and z. You're done!





