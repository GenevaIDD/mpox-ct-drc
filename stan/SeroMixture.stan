data {

 int N;          // N individuals
 vector[N] y;    // titer data
 int predL; // length of titer fit values
 vector[predL] y_fit; // titer fit values
 
}


parameters {

 real<lower=0, upper=1> sero;
 real mu0;              // mean seroneg
 real <lower=0> mu1;   // mean seropos
 real <lower=0> sd0;   // sd seroneg
 real <lower=0> sd1;   // sd seropos

}


transformed parameters {

  array[2] vector[N] pC;
  vector[N] log_lik;


  //--- likelihood calculation ---//
  for(n in 1:N){
    pC[1,n] = log(1 - sero) + normal_lpdf(y[n] | mu0, sd0);
    pC[2,n] = log(sero) + normal_lpdf(y[n] | mu0+mu1, sd1);
    log_lik[n] = log_sum_exp(pC[,n]);
  }
}


model {

 // priors
 sero ~ beta(1,1.5);
 mu0 ~ normal(4, 1);
 mu1 ~ normal(6, 1);
 sd0 ~ normal(0.5, 0.5);
 sd1 ~ normal(0.5, 0.5);

 // log-likelihood
 target += sum(log_lik);

}


generated quantities {

  // posterior predictive draws and component membership for each individual
  array[N] int  z;      // latent component: 0 = seroneg, 1 = seropos
  vector[N] y_rep;  // posterior predictive titer
  vector[predL] fitNeg;
  vector[predL] fitPos;
  vector[predL] fitAll;


  for(n in 1:N){
    
    real prob1 = exp(pC[2,n] - log_lik[n]);


    z[n] = binomial_rng(1, prob1);
    y_rep[n] = z[n] == 0 ? normal_rng(mu0, sd0)
                         : normal_rng(mu0 + mu1, sd1);
  }
  
  for(i in 1:predL){
    fitNeg[i] = (1 - sero) * exp(normal_lpdf(y_fit[i] | mu0, sd0));
    fitPos[i] = sero       * exp(normal_lpdf(y_fit[i] | mu0+mu1, sd1));
    fitAll[i] = fitNeg[i] + fitPos[i];
  }

}

