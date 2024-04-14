// Model for a single region, all years, no connection between years
//fitting as density for the observed number of deaths for different ages
//adding cohort effect for specific cohorts, common cohort effects across all years

functions {
  //return an int - to be used as index - from a real
  //https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/5
  int bin_search(real x, int min_val, int max_val){
    // This assumes that min_val >= 0 is the minimum integer in range, 
    //  max_val > min_val,
    // and that x has already been rounded. 
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1) %/% 2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1) %/% 2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }
  
  int[] signum_vec(vector x){ //sign function
    int n = num_elements(x);
    int res[n];
    for(i in 1:n){
      res[i] = x[i]<0 ? -1 : 1;
    }
    return res;
  }
  
  int signum(real x){ //sign function
    return x<0 ? -1 : 1;
  }
  
   real c_sn(real gamma){ //calculates c for skew normal reparametrization from centered to direct
    return signum(gamma) * cbrt(2*abs(gamma) / (4-pi()));
  }

  real[] c_sn_vec(real[] gamma){ //calculates c for skew normal reparametrization from centered to direct
    int n = num_elements(gamma);
    real res[n];
    for (i in 1:n){
      res[i]=signum(gamma[i]) * cbrt(2*abs(gamma[i]) / (4-pi()));
    }
    return res;
  }

  real[] mu_z_sn_vec(real[] c){ //calculates \mu_z for skew normal reparametrization from centered to direct
    int n = num_elements(c);
    real res[n];
    for (i in 1:n){
      res[i]=c[i]/sqrt(1+pow(c[i], 2));
    }
    return res;
  }
  
  real mu_z_sn(real c){ //calculates \mu_z for skew normal reparametrization from centered to direct
    return c/sqrt(1+pow(c, 2));
  }
  
  real[] xi_sn_vec(real[] mu, real[] sigma, real[] mu_z){ //calculates \mu_z for skew normal reparametrization from centered to direct
    int n = num_elements(mu);
    real res[n];
    for (i in 1:n){
      res[i]=mu[i]-sigma[i]*mu_z[i]/sqrt(1-pow(mu_z[i], 2));
    }
    return res;
  }
  
  real xi_sn(real mu, real sigma, real mu_z){
    return mu-sigma*mu_z/sqrt(1-pow(mu_z, 2));
  }
  
  real[] omega_sn_vec(real[] sigma, real[] mu_z){ 
    int n = num_elements(sigma);
    real res[n];
    for (i in 1:n){
      res[i]=sigma[i]/sqrt(1-pow(mu_z[i], 2));
    }
    return res;
  }
  
  real omega_sn(real sigma, real mu_z){
    return sigma/sqrt(1-pow(mu_z, 2));
  }

  real[] lambda_sn_vec(real[] mu_z){ 
    int n = num_elements(mu_z);
    real res[n];
    for (i in 1:n){
      res[i]=mu_z[i]*sqrt(pi()/2)/sqrt(1-pi()*pow(mu_z[i], 2)/2);
    }
    return res;
  }

  real lambda_sn(real mu_z){
    return mu_z*sqrt(pi()/2)/sqrt(1-pi()*pow(mu_z, 2)/2);
  }
}

// The input data is a matrix of life table deaths 'dx' of sum 'N'.
data {
  int<lower=0> x1; //minimum age
  int<lower=0> x2; //maximum age
  int<lower=0> t1; //starting year
  int<lower=0> t2; //final year
  int<lower=0> N[t2-t1+1]; //population size
  int<lower=0> dx[t2-t1+1, x2-x1+1]; //life table deaths
  int<lower=0> N_coh; //number of cohort effects
  int<lower=0> cohorts[N_coh]; //cohort year of birth
}

transformed data {
  int ages[x2-x1+1];
  for (y in 1:(x2-x1+1)){
    ages[y]=x1+y-1;
  }
  
  real ages2[x2-x1+1]; //mid-year (for density-like estimation)
  int coh_indicator[t2-t1+1, x2-x1+1]; //indicator of cohort effect - gives either 0 or index of cohort effect
  for (y in 1:(x2-x1+1)){
    if(x1==0 && y==1){
      ages2[y]=0.1; //deaths at age 0 are closer to birth on average
    } else {
      ages2[y]=x1+y-0.5; //deaths at other ages can be considered uniformly distributed
    }
  }
  for (t in 1:(t2-t1+1)){
    for (y in 1:(x2-x1+1)){
      coh_indicator[t, y]=0;
    }
    for (c in 1:N_coh)
    {
      //print("t: ", t, ", c: ", c, ", cohorts[c]: ", cohorts[c], ", t-cohorts[c]: ", t-cohorts[c]);
      if(t1+t-1-cohorts[c]>=x1 && t1+t-1-cohorts[c]<=x2){ //cohort effect relevant for given age range
      coh_indicator[t, t1+t-cohorts[c]]=c;
      }
    }
  }
}

/* See Zanotto et al., A Mixture‑Function Mortality Model:
Illustration of the Evolution of Premature Mortality (2021)
Three components
*/
parameters {
  simplex[3] zeta[t2-t1+1]; //3 mixing parameters (infant, premature, adult mortality)

  //infant mortality (1 parameter)
  real<lower=0> sigma_i[t2-t1+1]; //scale
  
  //premature mortality (3 parameters)

  positive_ordered[2] mu[t2-t1+1]; //location for premature and adult mortality
  real<lower=0> sigma_m[t2-t1+1]; //scale
  real <lower=0, upper=0.99527> gamma_m[t2-t1+1]; //shape 

  //adult mortality (3 parameters)
  real<lower=0> sigma_M[t2-t1+1]; //scale
  real <lower=-0.99527, upper=0.99527> gamma_M[t2-t1+1]; //shape

  real log_alpha_raw[N_coh]; //cohort effect on log-density
  real<lower=0, upper=1> alpha_mu; //mean effect for alpha_raw
}

transformed parameters {
  
  //direct reparametrization
  real xi_m[t2-t1+1];
  real<lower=0> omega_m[t2-t1+1];
  real lambda_m[t2-t1+1];
  
  real xi_M[t2-t1+1];
  real<lower=0> omega_M[t2-t1+1];
  real lambda_M[t2-t1+1];
  
  real denominator_m[t2-t1+1];
  real denominator_M[t2-t1+1];
  
  real denominator[t2-t1+1];
  
  real log_p_x[t2-t1+1, x2-x1+1]; //truncated mixture
  real log_alpha[N_coh];
  real log_alpha_mean = mean(log_alpha_raw);
  
  for (c in 1:N_coh){
    log_alpha[c]=log_alpha_raw[c]-log_alpha_mean;
  }
  
  for (t in 1:(t2-t1+1)){
    xi_m[t] = xi_sn(mu[t,1], sigma_m[t], mu_z_sn(c_sn(gamma_m[t]))); //location 
    omega_m[t] = omega_sn(sigma_m[t], mu_z_sn(c_sn(gamma_m[t]))); //scale
    lambda_m[t] = lambda_sn(mu_z_sn(c_sn(gamma_m[t]))); //shape

    xi_M[t] = xi_sn(mu[t,2], sigma_M[t], mu_z_sn(c_sn(gamma_M[t]))); //location
    omega_M[t] = omega_sn(sigma_M[t], mu_z_sn(c_sn(gamma_M[t]))); //scale
    lambda_M[t] = lambda_sn(mu_z_sn(c_sn(gamma_M[t]))); //shape
    
    denominator_m[t]=(skew_normal_cdf(ages[x2-x1 +1]+1, xi_m[t], omega_m[t], lambda_m[t])-skew_normal_cdf(ages[1], xi_m[t], omega_m[t], lambda_m[t]));
    denominator_M[t]=(skew_normal_cdf(ages[x2-x1 +1]+1, xi_M[t], omega_M[t], lambda_M[t])-skew_normal_cdf(ages[1], xi_M[t], omega_M[t], lambda_M[t]));
    
    denominator[t]=negative_infinity();

    for(x in 1:(x2-x1 +1)){
      log_p_x[t, x]=log_sum_exp( //f(x;\theta)
      log_sum_exp(
        log(zeta[t, 1]) + log(2) + normal_lpdf(ages2[x] | 0, sigma_i[t]), //log(2) added since it's half-normal and the probability is doubled
        log(zeta[t, 2]) + skew_normal_lpdf(ages2[x] | xi_m[t], omega_m[t], lambda_m[t]) - log(denominator_m[t])),
        log(zeta[t, 3]) + skew_normal_lpdf(ages2[x] | xi_M[t], omega_M[t], lambda_M[t]) - log(denominator_M[t])
        );
      if(coh_indicator[t, x]>0){
        denominator[t]=log_sum_exp(denominator[t], log_alpha[coh_indicator[t, x]] + log_p_x[t, x]);
      } else {
        denominator[t]=log_sum_exp(denominator[t], log_p_x[t, x]); //year without a cohort effect (hence alpha fixed at 1, hence log=0)
      }
    }
    //print("log_alpha: ", log_alpha, ", p_x_t: ", log_p_x[t], ", denominator: ", denominator[t]);
  }
}

// The model to be estimated.
model {
  alpha_mu~beta_proportion(0.5, 3); //slightly informative, slight shrinking towards 0
  
  for(c in 1:N_coh){
    log_alpha_raw[c] ~ uniform(log(alpha_mu+0.5)-0.7, log(alpha_mu+0.5)+0.4); //very wide for realistic cohort effects (50%, 150%)
  }
  
  for (t in 1:(t2-t1+1)){
    //mixing parameter
    zeta[t] ~ dirichlet([1, 1, 8]);
    
    //premature mortality
    mu[t,1] ~ normal(30, 5); 
    sigma_m[t] ~ inv_gamma(0.001, 0.001); 
    gamma_m[t] ~ normal(0, 1); 
    
    //adult mortality
    mu[t,2] ~ normal(80, 10);
    sigma_M[t] ~ inv_gamma(0.001, 0.001); 
    gamma_M[t] ~ normal(0, 1);
    
    
    for(x in 1:(x2-x1 +1)){
      if(coh_indicator[t, x]>0){
        target+=dx[t, x]*(
          log_alpha[coh_indicator[t, x]] + //alpha_{c}
          log_p_x[t, x] - denominator[t]);
      } else {
        target+=dx[t, x]*(
          log_p_x[t, x] - denominator[t]);
      }
    }
  }

}

generated quantities {
 
  simplex[x2-x1+1] theta[t2-t1+1]; //probabilities of death
  
  for(t in 1:(t2-t1 +1)){
    for(x in 1:(x2-x1 +1)){
      if(coh_indicator[t, x]){
        theta[t,x]=exp(log_alpha[coh_indicator[t, x]] +
        log_p_x[t, x] - denominator[t]
        );
      } else {
        theta[t,x]=exp(log_p_x[t, x] - denominator[t]
        );
      }
    }
    //print("theta t: ", theta[t], ", sum of theta t: ", sum(theta[t]));
    theta[t] = theta[t]/sum(theta[t]); //this should not be necessary
  }

  real<lower=0> dx_rep[t2-t1+1, x2-x1+1];
  for(t in 1:(t2-t1+1)){
    for(x in 1:(x2-x1+1)){
      dx_rep[t, x] = theta[t, x]*N[t];
    }
  }

}
