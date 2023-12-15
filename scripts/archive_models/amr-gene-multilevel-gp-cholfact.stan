// random intercept per participant
// Time-varying effect of drug effect  
// binary response variable

functions {
  //generate covariance matrix
  matrix generate_K(matrix dist_sq_matrix, int n_samples, real alpha_sq, real sigma_sq, real length_scale_sq) {
    matrix[n_samples, n_samples] out_matrix;
 //   print("gen fn");
    for (r in 1:n_samples) {
      for (c in r:n_samples) {
//        print("r,c,val");
//        print(r);
//        print(c);
//        print(dist_sq_matrix[r,c]);
        if (dist_sq_matrix[r,c] < 0) {
          if (r==c) {
            out_matrix[r,c] = alpha_sq * exp(-dist_sq_matrix[r,c] / (2 * length_scale_sq)) + sigma_sq + 1e-9;
          }
          else {
            out_matrix[r,c] = 0; 
            out_matrix[c,r] = out_matrix[r,c];
          }
        }
        else {
          if (r==c) {
            out_matrix[r,c] = alpha_sq * exp(-dist_sq_matrix[r,c] / (2 * length_scale_sq)) + sigma_sq + 1e-9;
          }
          else {
            out_matrix[r,c] = alpha_sq * exp(-dist_sq_matrix[r,c] / (2 * length_scale_sq));
            out_matrix[c,r] = out_matrix[r,c];
          }
        }
 //       print("out");
 //       print(out_matrix[r,c]);
      }
    }
    
 // print(out_matrix);
  return out_matrix;
  }



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
array[n] int<lower=0, upper=1> y;                    // response variable [0,1]
matrix[n,n] dist_matrix;         // distance matrix encoding correlations
}

parameters{
//-----------participant effects ---------------------------------------------
real<lower=0> sigma;       // sd of correlation at t = 0
real<lower=0> alpha;
real<lower=0> length_scale;         // length scale for quadratic kernal
// ----------drug effects ----------------------------------------
vector[n_covariates] beta;    // effect on gene presence of drug exposure
real<lower=0> tau;            // mean lifetime of drug effect- same for all drugs
// ------------- intercept ------------------------------------------
real beta_0;        // intercept
vector[n] eta;
}

transformed parameters{
}

model{
// generate covaraiance matrix
matrix[n,n] K;
matrix[n,n] L_K;
vector[n] f;
K = generate_K(dist_matrix, n, alpha^2, sigma^2, length_scale^2);
//print("K:");
//print(K);
L_K = cholesky_decompose(K);
f = L_K * eta;
//print("gamma");
//print(gamma);

// likelihood -----------------------------------------------------------------
// iterate over samples
for (j in 1:n) {
  y[j] ~ bernoulli_logit(beta_0 + f[j] + (effect_vector(n_covariates, t_e[j], tau) * beta));
}
// normally distributed individual effects ------------------------------------
//priors - weakly informative -------------------------------------------------
beta ~ student_t(3, 0, 3);
beta_0 ~ student_t(3, 0, 3);
tau ~ student_t(3,0,50);
sigma ~ normal(0,3);
alpha ~ normal(0,3);
length_scale ~ inv_gamma(2,2);
eta ~ std_normal();
}

generated quantities {
}

