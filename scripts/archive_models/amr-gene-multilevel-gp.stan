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
real<lower=0> sigma_sq;       // sd of correlation at t = 0
real<lower=0> alpha_sq;
real<lower=0> length_scale_sq;         // length scale for quadratic kernal
// ----------drug effects ----------------------------------------
vector[n_covariates] beta;    // effect on gene presence of drug exposure
real<lower=0> tau;            // mean lifetime of drug effect- same for all drugs
// ------------- intercept ------------------------------------------
real beta_0;        // intercept
vector[n] gamma;    // gaussian process bit
}

transformed parameters{
}

model{
// generate covaraiance matrix
matrix[n,n] K;

// matrix[n,n] L_K;
K = generate_K(dist_matrix, n, alpha_sq, sigma_sq, length_scale_sq);
//print("K:");
//print(K);
L_K = cholesky_decompose(K);
gamma ~ multi_normal_cholesky(rep_vector(0,n), K);
//print("gamma");
//print(gamma);

// likelihood -----------------------------------------------------------------
// iterate over samples
for (j in 1:n) {
  y[j] ~ bernoulli_logit(beta_0 + gamma[j] + (effect_vector(n_covariates, t_e[j], tau) * beta));
}
// normally distributed individual effects ------------------------------------
//priors - weakly informative -------------------------------------------------
beta ~ student_t(3, 0, 3);
beta_0 ~ student_t(3, 0, 3);
tau ~ student_t(3,0,50);
sigma_sq ~ normal(0,3);
alpha_sq ~ normal(0,3);
length_scale_sq ~ gamma(4, 0.125);
}

generated quantities {
}

