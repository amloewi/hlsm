
data {
    int<lower=0> N;   // number of nodes
    int<lower=0> edges[N,N]; //for now non-blank => 1
    real<lower=0> sigma_alpha; // alpha is the density
    vector[2] mu_z;             // mean of the latent positions
    matrix<lower=0>[2,2] sigma_z; // variance of the latent positions
    matrix<lower=0, upper=1>[2,2] sigma_fixed;
}

parameters {
    real alpha;     // density intercept
    vector[2] z[N]; // latent positions
}


model {

    real d; // just a convenience variable, to split up a long line

    // prior on alpha
    alpha ~ normal(0, sigma_alpha);

    // latent variables
    z[1] ~ multi_normal(mu_z, sigma_fixed);
    for(i in 2:N) {
        z[i] ~ multi_normal(mu_z, sigma_z);
    }

    for(i in 1:N){
        for(j in 1:N){
            // This works with euclidean^2, NOT with sqrt
            # d = (pow(z[i,1]-z[j,1],2) + pow(z[i,2]-z[j,2],2));
            d = dot_self(z[i] - z[j])
            edges[i,j] ~ bernoulli_logit(alpha - d);
        }
    }
}
