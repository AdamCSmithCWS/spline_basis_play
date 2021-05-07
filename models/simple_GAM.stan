
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  int<lower=2> nknots;
  matrix[N,nknots] Xb;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //vector[N] mu;
  real alpha;
  real<lower=0> sd_beta;
  real<lower=0> sigma;
  vector[nknots] betas_raw;
}

transformed parameters {
  
  vector[nknots] betas;
  vector[N] mu;
  
  
betas = sd_beta*betas_raw;
 mu = Xb*betas + alpha;

}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  sigma ~ student_t(3,0,2.5);
  sd_beta ~ student_t(3,0,2.5);
  
  alpha ~ student_t(3,0,2.5);
  betas_raw ~ normal(0,1);
  sum(betas_raw) ~ normal(0,0.001*nknots);
  
  
  x ~ normal(mu, sigma);
}

