// random intercept per participant
// Time-varying effect of drug effect  
// binary response variable

functions {
  // return vector of covariate effect to dot product with betas
  row_vector effect_vector(int vector_length, row_vector t_e, real tau) {
  row_vector[vector_length] out_row_vector;
  for (i in 1:vector_length) {  
    if(t_e[i] < 0) {
        out_row_vector[i] = 0;
    }
    else {
      out_row_vector[i] = exp(-t_e[i] / tau);
    }
  }
    return out_row_vector;
  }
}

data {
int<lower=1> n;                 //total number of observations
int n_participants;             // number of participants
int n_covariates;               // number of covariates
array[n] int pid;               // map observation on to participant id
array[n] row_vector[n_covariates] t_e;                  // days since drug exposure; <0 = no exposure
array[n] int<lower=0> y;                    // response variable - int
array[n] real<lower=0> d;                     // log sequencing depth for obs n
}


parameters{
//-----------participant effects ---------------------------------------------
vector<lower=0>[n_participants] gamma;    // intercept for each participant
real<lower=0> sigma_participant;       // sd of participant intercepts
// ----------drug effects ----------------------------------------
vector[n_covariates] beta;    // effect of drug exposure
real<lower=0> tau;            // mean lifetime of drug effect- same for all drugs
// ------------- intercept ------------------------------------------
real<lower=0> alpha;
// distribution effects
real<lower=0> phi;         // overdispertion parameter;
}

transformed parameters{
}

model{
// likelihood -----------------------------------------------------------------
for (j in 1:n) {
  y[j] ~ neg_binomial_2_log(alpha + gamma[pid[j]] + d[j] + (effect_vector(n_covariates,
  t_e[j], tau) * beta), phi);
}
// normally distributed individual effects ------------------------------------
gamma ~ normal(0, sigma_participant);
//priors - weakly informative -------------------------------------------------
beta ~ student_t(3, 0, 3);
alpha ~ student_t(3, -3, 3);
sigma_participant ~ student_t(3,0,3);
tau ~ student_t(3,0,50);
phi ~ student_t(3, 0, 3);
}

generated quantities {

array[n] int<lower=0> y_pred;

for (j in 1:n) {
  y_pred[j] = neg_binomial_2_log_rng(alpha + gamma[pid[j]] + d[j] + (effect_vector(n_covariates,
  t_e[j], tau) * beta), phi);
}
}

