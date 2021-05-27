// This is a Stan implementation of the bbsBayes gamye model
// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues




data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations        
  int<lower=1> year[ncounts]; // year index
  int<lower=1> route[ncounts]; // route index
 
 // real g1; //first gamma prior parameter
 // real g2; //first gamma prior parameter
 

  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basispred; // basis function matrix
 real sdB_pr; //prior on sd of t-dist prior for sdBETA

}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
 
 
  vector[nyears] yeareffect_raw;

  vector[nroutes] rte_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdrte;    // sd of route effects
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=0> sdyear;    // sd of year effects
 // real<lower=4,upper=500> nu; // df of t-distribution > 4 so that it has a finite mean, variance, kurtosis
  
  vector[nknots_year] BETA_raw;//_raw; 
 real alpha; // overall intercept
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] rte; // route intercepts
  vector[nyears] Y_pred;  

  vector[nyears] yeareffect;
  vector[ncounts] noise;             // over-dispersion
  vector[nknots_year] BETA;
  
  
  BETA = sdBETA*BETA_raw;
  
  Y_pred = year_basispred * BETA; 
  
 
     rte = sdrte*rte_raw;

yeareffect = sdyear*yeareffect_raw;
// intercepts and slopes

  noise = sdnoise*noise_raw;
  
  

  for(i in 1:ncounts){
    E[i] =  alpha + Y_pred[year[i]] + yeareffect[year[i]] + rte[route[i]] + noise[i];
  }
  
  }
  
  
  
model {
  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1);// ~ student_t(4,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
   
  sdrte ~ normal(0,1); //prior on sd of route effects
  sdyear ~ gamma(2,2); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
  //sdBETA ~ gamma(2,2); // prior on sd of GAM parameters
  sdBETA ~ normal(0,sdB_pr); // prior on sd of GAM parameters
  //sdBETA ~ student_t(3,0,sdB_pr); // prior on sd of GAM parameters
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution


  rte_raw ~ std_normal();//route effects
  sum(rte_raw) ~ normal(0,0.001*nroutes);
 
 
  yeareffect_raw ~ std_normal();
  //soft sum to zero constraint on year effects within a stratum
  sum(yeareffect_raw) ~ normal(0,0.001*nyears);
  
 alpha ~ std_normal();
  
  BETA_raw ~ std_normal();// prior on fixed effect mean GAM parameters
  //sum to zero constraint
  //sum(BETA_raw) ~ normal(0,0.001*nknots_year);
  
  

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation


}

 generated quantities {

   real<lower=0> n[nyears];
   real<lower=0> nsmooth[nyears];
   real<lower=0> retrans_noise;
  vector[ncounts] log_lik;
  real retrans_yr = 0.5*(sdyear^2); //retransformation factor for sdyear

  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  
retrans_noise = 0.5*(sdnoise^2);//retransformation factor for sdnoise

for(y in 1:nyears){

 
  real n_t[nroutes];
  real nsmooth_t[nroutes];

        for(t in 1:nroutes){
            
      n_t[t] = exp(alpha + Y_pred[y] + rte[t] + yeareffect[y] + retrans_noise);
      nsmooth_t[t] = exp(alpha + Y_pred[y] + rte[t] + retrans_noise + retrans_yr);
        }
        n[y] = mean(n_t);
        nsmooth[y] = mean(nsmooth_t);



    }
  }





