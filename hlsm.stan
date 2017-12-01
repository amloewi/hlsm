
data {
    int<lower=0> N;   // number of nodes
    int<lower=0> K;   // number of layers in the multigraph
    int<lower=0, upper=1> edges[N,N,K]; //for now non-blank => 1
    real<lower=0> sigma_alpha; // alpha is the layer intercept
    vector[2] mu_b; // prior mean for the base positions
    matrix<lower=0>[2,2] sigma_b;     // variance of the base positions -- number?
    matrix<lower=0>[2,2] sigma_z[K];     // variance of the layer positions FROM b;
                                  // should be K so that you can penalize LAYERS
}

parameters {
    // Must be some variant of vector ~ mvn(vector, matrix)

    real alpha[K];    // layer intercepts; scalar
    // Strange syntax; vector[2] is what the structure holds,
    // x[A,B] are the dimensions of the structure
    vector[2] b[N];   // base positions; in 2d
    vector[2] z[N,K]; // layer positions; [person, layer, x/y coord]
    // DO I NEED A (fourth now) DIMENSION TO THIS? mvlsm has one.
}


model {

    real d; // just a convenience variable, to split up a long line

    // prior on alpha
    for (k in 1:K) {
        alpha[k] ~ normal(0, sigma_alpha);
    }

    // prior on the base positions
    for(i in 1:N){
        b[i] ~ multi_normal(mu_b, sigma_b); # XXX
    }

    // priors on the layer positions
    for (i in 1:N) {
        for (k in 1:K) {
            z[i,k] ~ multi_normal(b[i], sigma_z[k]);
            // a laplacian prior is like a lasso on the variances, meaning
            // this could zero out a whole LAYER (model a layer as simply
            // being equal to the base positions)
            // z[i,k] ~ double_exponential(b[i], sigma_z_vector[k]);
        }
    }

    # Fitting the parameters to the edges
    for(i in 1:N){
        for(j in 1:N){
            for(k in 1:K){
                d = dot_self(z[i,k] - z[j,k]);
                edges[i,j,k] ~ bernoulli_logit(alpha[k] - d);

                #    y[n] ~ categorical(softmax(beta * x[n]));

            }
        }
    }
}
