//----- Null (single-component) model -----//
// All Ct>0 observations are modelled as true infections: y ~ Normal(mu, sdI).
// No latent class mixture, no environmental (false-positive) component.

data {

  int<lower=0> N;                    // N Ct results from cases
  vector[N] y;                       // log Ct data from cases
  vector[N] t;                       // time since symp onset
  array[N] int site;                 // location
  vector[4] NperSite;                // N per location
  vector[45] ct_range;               // range for posthoc summaries
  int<lower=1> NCov;                 // N covariates
  array[NCov] int NG;                // N groups per covariate
  array[N,NCov] int cov;             // covariate data
  array[NCov,2] int ind_cov;         // covariate indices
  int<lower=1> K;                    // N spline basis functions
  matrix[N,K] B;                     // B-spline basis matrix

}


parameters {

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
  vector[N] log_lik;                 // per-observation log-likelihood
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

  for(n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | mu[n], sdI);

  }
}


model {

  // priors
  c0 ~ normal(3,1);                  // weakly informative prior on log Ct intercept
  sdI ~ normal(0.2,0.05);            // prior on infection Ct sd
  alpha_raw ~ normal(0,1);           // weakly regularising prior on covariate group effects
  tau ~ exponential(1);              // shrinkage prior on per-covariate effect scale
  sigma_beta ~ normal(1,0.1);        // smoothing penalty scale

  beta[1] ~ normal(0, 5);            // flat priors on first two spline coefficients
  if(K >= 2) beta[2] ~ normal(0, 5);
  if(K >= 3)
    for(k in 3:K)
      (beta[k] - 2*beta[k-1] + beta[k-2]) ~ normal(0, sigma_beta); // second-difference smoothing penalty

  y ~ normal(mu, sdI);

}


generated quantities {

  vector[45] pI;                     // marginal infection density over ct_range
  matrix[NCov,max(NG)] alpha_eff = rep_matrix(0, NCov, max(NG)); // scaled covariate effects: alpha * tau

  // pI: infection density evaluated on ct_range grid; untruncated to match the likelihood
  real mean_mu = mean(mu);
  for(i in 1:45)
    pI[i] = exp(normal_lpdf(ct_range[i] | mean_mu, sdI));

  for(c in 1:NCov)
    alpha_eff[c,1:NG[c]] = alpha[c,1:NG[c]] .* tau[c];

  // y_rep: draws from Normal(mu[n], sdI) — no truncation or mixture needed
  vector[N] y_rep;
  for(n in 1:N)
    y_rep[n] = normal_rng(mu[n], sdI);

}
