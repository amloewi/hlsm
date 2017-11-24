
data {
    int<lower=0> N;   // number of nodes
    int<lower=0> K;   // number of layers in the multigraph
    int<lower=0, upper=1> edges[N,N,K]; //for now non-blank => 1
    real<lower=0> sigma_alpha; // alpha is the layer intercept
    vector[2] mu_b;
    matrix<lower=0>[2,2] sigma_b;     // variance of the base positions -- number?
    matrix<lower=0>[2,2] sigma_z[K];     // variance of the layer positions FROM b;
                                  // should be K so that you can penalize LAYERS
}

parameters {
    // Must be some combination of vector ~ mvn(vector, matrix)

    real alpha[K];    // layer intercepts; scalar
    // XXX NOT SURE IF THIS IS WHAT I WANT; vector[2]?
    vector[2] b[N];   // base positions; in 2d
    vector[2] z[N,K]; // layer positions; person, layer, x/y coord
    // DO I NEED A (fourth now) DIMENSION TO THIS? mvlsm has one.
}


model {

    real d; // just a convenience variable, to split up a long line

    // prior on alpha
    for (k in 1:K) {
        alpha[k] ~ normal(0, sigma_alpha);
    }

    // prior on the base positions -- IS IT? WILL THIS FREEZE IT?
    for(i in 1:N){
        b[i] ~ multi_normal(mu_b, sigma_b); # XXX
    }

    for (i in 1:N) {
        for (k in 1:K) {
            z[i,k] ~ multi_normal(b[i], sigma_z[k]);
            // z[i,k] ~ double_exponential(b[i], sigma_z[k]); // requires VEC sig
        }
    }

    for(i in 1:N){
        for(j in 1:N){
            for(k in 1:K){
                // pow required a real and this is easy, if messy
                d = (pow(z[i,k,1]-z[j,k,1],2) + pow(z[i,k,2]-z[j,k,2],2));
                edges[i,j,k] ~ bernoulli_logit(alpha[k] - d);
            }
        }
    }
}
