soft_thresh = function(x, g) {
    x = as.vector(x)
    w1 = which(g >= abs(x))
    w2 = which(g < abs(x) & x > 0)
    w3 = which(g < abs(x) & x < 0)
    ret = x
    ret[w1] = 0
    ret[w2] = x[w2]-g
    ret[w3] = x[w3]+g
    ret
}
glmnet_ref = function(X, y, lambda, alpha, family=binomial, maxit=10, tol=1e-08)
{
    beta = matrix(rep(0,ncol(X)), ncol=1)
    for(j in 1:maxit)
    {
        beta_outer_old = beta
        eta    = as.matrix(X %*% beta)
        g      = family()$linkinv(eta)
        gprime = family()$mu.eta(eta)
        z      = eta + (y - g) / gprime
        W      = as.vector(gprime^2 / family()$variance(g))
        wx_norm = colSums(W*X^2)    
        for (k in 1:maxit) {
            beta_inner_old = beta
            for (l in 1:length(beta)) {
                beta[l] = soft_thresh(sum(W*X[,l]*(z - X[,-l] %*% beta_inner_old[-l])),
                                      sum(W)*lambda*alpha)
            }
            beta = beta / (wx_norm + lambda*(1-alpha))
            if(sqrt(as.double(crossprod(beta-beta_inner_old))) < tol) break    }
        if (sqrt(as.double(crossprod(beta-beta_outer_old))) < tol) break
    }
    list(beta=beta,iterations=j)
}

# compare with glmnet to test
library(glmnet)
library(datasets)
library(mlbench)
data(Sonar)

X <- as.matrix(Sonar[,-ncol(Sonar)])
y <- Sonar$Class
sonar.glm <- glmnet(X, y=="R", family="binomial")
lambda <- sonar.glm$lambda[10]
coef(sonar.glm, s=lambda)

sonar.hack <- glmnet_ref(X, y=="R", lambda, 1, maxit=20) # so -- it doesn't work.
sonar.hack$beta
