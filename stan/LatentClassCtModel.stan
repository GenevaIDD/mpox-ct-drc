//----- Latent class model with covariate effects -----//


data {

  int<lower=0> N;                    // N Ct results from cases
  int<lower=0> NE;                   // N environmental samples
  int<lower=1> N_sites;              // number of sites
  vector[N] y;                       // log Ct data from cases
  vector[N] t;                       // time since symp onset
  vector[NE] yE;                     // log environmental Ct data
  array[N] int site;                 // location (integer 1..N_sites)
  vector[N_sites] NperSite;          // N per location
  vector[N_sites] N_ct_thresh;       // N samples with Ct < ct_threshold per site (for truepI_all weighting)
  real<lower=0> ct_threshold;        // Ct positivity threshold (raw scale)
  vector[45] ct_range;               // range for posthoc summaries
  int<lower=1> NCov;                 // N covariates
  array[NCov] int NG;                // N groups per covariate
  array[N,NCov] int cov;             // covariate data
  array[NCov,2] int ind_cov;         // covariate indices
  int<lower=1> K;                    // N spline basis functions
  matrix[N,K] B;                     // B-spline basis matrix

}


parameters {

  vector<lower=0, upper=1>[N_sites] piI;  // prop of true infections among Ct data by site
  real muE;                          // mean of environmental Ct component
  real<lower=0> sdE;                 // sd of environmental Ct component
  real<lower=0> sdI;                 // sd of true infection Ct component
  real c0;                           // intercept
  vector[K] beta;                    // spline coefficients
  real<lower=0> sigma_beta;          // smoothing penalty
  vector[sum(NG)] alpha_raw;         // uncentered group effects of each covariate level
  vector<lower=0>[NCov] tau;         // per-covariate effect scale

}


transformed parameters {

  vector[N] mu;                      // expected mean log Ct per observation
  matrix[NCov,max(NG)] alpha = rep_matrix(0, NCov, max(NG)); // centred covariate group effects
  array[N,2] real pC;                // log component densities [1]=infection, [2]=environmental
  vector[N] log_lik;                 // per-observation log-likelihood
  real sumloglik;                    // sum of log-likelihoods (added to target)
  {
  int ix = 1;
  vector[N] delta = rep_vector(0, N); // total covariate offset on log Ct
  for(c in 1:NCov) for(g in 1:NG[c]){
    alpha[c,g] = alpha_raw[ix] - mean(alpha_raw[ind_cov[c,1]:ind_cov[c,2]]); // centre within covariate
    ix = ix + 1;
  }

  for(n in 1:N)
    for(c in 1:NCov) delta[n] += alpha[c,cov[n,c]] * tau[c];

  mu = c0 + B*beta + delta;  // intercept + spline trajectory + covariate offset

  real log_thresh = log(ct_threshold);  // truncation point

  for(n in 1:N){
    // Infection component: Normal(mu[n], sdI) truncated at log(40) on the right.
    // For y > log(ct_threshold): log-density = -Inf (zero probability of being a true infection).
    // For y <= log(ct_threshold): subtract log-normalising constant log(Phi((log_thresh-mu)/sdI)).
    real log_p_inf;
    if(y[n] <= log_thresh)
      log_p_inf = normal_lpdf(y[n] | mu[n], sdI) - normal_lcdf(log_thresh | mu[n], sdI);
    else
      log_p_inf = negative_infinity();

    pC[n,1] = log(piI[site[n]])   + log_p_inf;
    pC[n,2] = log1m(piI[site[n]]) + normal_lpdf(y[n] | muE, sdE);
    log_lik[n] = log_sum_exp(pC[n,]);
  }
  sumloglik = sum(log_lik);

  }
}


model {

  // priors
  muE ~ normal(mean(yE), 0.003);     // anchored to sample mean of environmental Ct
  sdE ~ normal(sd(yE), 0.003);       // anchored to sample sd of environmental Ct
  piI ~ beta(10,1);                 
  c0 ~ normal(3,1);                 
  sdI ~ normal(0.2,0.05);            
  alpha_raw ~ normal(0,1);           
  tau ~ exponential(1);              
  sigma_beta ~ normal(1,0.1);        

  beta[1] ~ normal(0, 5);           
  if(K >= 2) beta[2] ~ normal(0, 5);
  if(K >= 3)
    for(k in 3:K)
      (beta[k] - 2*beta[k-1] + beta[k-2]) ~ normal(0, sigma_beta); // second-difference smoothing penalty

  target += sumloglik;

}


generated quantities {

  vector[N_sites] truepI;            // P(true infection | Ct < ct_threshold) by site
  matrix[NCov,max(NG)] alpha_eff = rep_matrix(0, NCov, max(NG)); // scaled covariate effects: alpha * tau
  real truepI_all;                   // P(true infection | Ct < ct_threshold) across all sites

  real log_thresh = log(ct_threshold);

  // truepI[s] = P(infection | Ct < ct_threshold, site s)
  for(s in 1:N_sites) truepI[s] = piI[s] / (N_ct_thresh[s]/NperSite[s]);
  truepI_all = sum(truepI .* N_ct_thresh) / sum(N_ct_thresh);

  // Mixture densities (pI, pEnv, pAll) are computed in R via get_marginal_dens(),
  // which correctly averages over individual mu values rather than using mean(mu).

  for(c in 1:NCov)
    alpha_eff[c,1:NG[c]] = alpha[c,1:NG[c]] .* tau[c];

  // y_rep: infection draws use the inverse-CDF method for a right-truncated
  // Normal: if U ~ Uniform(0, Phi((log_thresh - mu)/sdI)) then
  // mu + sdI * inv_Phi(U) ~ Normal(mu, sdI) truncated at log_thresh.
  vector[N] y_rep;
  for(n in 1:N){
    int z = bernoulli_rng(piI[site[n]]);
    if(z == 1){
      real ub_std = (log_thresh - mu[n]) / sdI;
      y_rep[n] = mu[n] + sdI * inv_Phi(uniform_rng(0.0, Phi(ub_std)));
    } else {
      y_rep[n] = normal_rng(muE, sdE);
    }
  }

}
