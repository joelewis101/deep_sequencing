// random intercept per participant
// Time-varying effect of drug effect  
// binary response variable

functions {

  // return matrix of covariate effect for product with betas
  matrix return_X(matrix X, real tau) {
    matrix[rows(X),cols(X)] out_X;
   // print("out_x");
    //print(out_X);
    for (r in 1:rows(X)) {
      for (c in 1:cols(X)) {
    //  print("r,c, X[r,c]");
    //  print(r);
    //  print(c);
    //  print(X[r,c]);
        if(X[r,c] < 0) {
          out_X[r,c] = 0;
        }
        else {
          out_X[r,c] = exp(-X[r,c] / tau);
        }
      }
    }
  //print(out_X);
  return(out_X);
  }
}

data {
int<lower=1> n;                 //total number of observations
array[n] real t;                      // day of observation
int n_participants;             // number of participants
array[n_participants] int N;           // number of obs per participant
int n_covariates;               // number of covariates
matrix[n,n_covariates] t_e;                  // days since drug exposure; <0 =
                                             //no exposure, one row per obs, one
                                             //col per drug
array[n] int<lower=0> y;                    // response variable - int > 0
array[n] real<lower=0> d;                     // log sequencing depth for obs n
}

parameters{
// ------ distribution parameters ----------------------------
real<lower=0> phi; // overdispertion parameter
//-----------participant effects ---------------------------------------------
//real<lower=0> sigma;       // sd of correlation at t = 0
real<lower=0> alpha;       // magnitude for quadratic kernal
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
// iterate over participants
  int p;  // pointer
  p = 1;

for (i in 1:n_participants) {
  matrix[N[i], N[i]] L_K; // cholesky decomposed covariance matrix for this participant
  vector[N[i]] f;        // GP effect for this participant
  matrix[N[i], n_covariates] X;   // covariate matrix for this participant
  vector[N[i]] d_p;  //depth for thi sparticipants samples
  //print("t");
  //print(t[p:(p+N[i]-1)]);
  L_K = cholesky_decompose( gp_exp_quad_cov(t[p:(p+N[i]-1)], alpha, length_scale));
  //L_K = cholesky_decompose( gp_exp_quad_cov(t[p:(p+N[i]-1)], alpha, length_scale) + diag_matrix(rep_vector(sigma, N[i])));
  //print("L_K");
  //print(L_K);
  f = L_K * eta[p:(p+N[i]-1)];
  //print("f");
  //print(f);
  //print("t_e");
  //print(t_e[p:(p+N[i]-1),]);
  X = return_X(t_e[p:(p+N[i]-1),], tau);
  //print("X");
  //print(X);
  d_p = to_vector(d[p:(p+N[i] -1)]);
  y[p:(p+N[i] -1)] ~ neg_binomial_2_log(beta_0 + f + X*beta + d_p, phi);
  p = p + N[i];
}


// normally distributed individual effects ------------------------------------
//priors - weakly informative -------------------------------------------------
beta ~ student_t(3, 0, 3);
beta_0 ~ student_t(3, -3, 3);
tau ~ student_t(3,0,3);
//sigma ~ student_t(3,0,3);
alpha ~ student_t(3,0,3);
//length_scale ~ inv_gamma(5,5);
length_scale  ~ inv_gamma(3,3);
eta ~ std_normal();
phi ~ student_t(3, 0, 3);
}

generated quantities {
  array[n] int<lower=0> y_sim;

  int p;  // pointer
  p = 1;

  for (i in 1:n_participants) {
    matrix[N[i], N[i]] L_K; // cholesky decomposed covariance matrix for this participant
    vector[N[i]] f;        // GP effect for this participant
    matrix[N[i], n_covariates] X;   // covariate matrix for this participant
    vector[N[i]] d_p;  //depth for thi sparticipants samples
    //print("t");
    //print(t[p:(p+N[i]-1)]);
    L_K = cholesky_decompose( gp_exp_quad_cov(t[p:(p+N[i]-1)], alpha, length_scale));
    //L_K = cholesky_decompose( gp_exp_quad_cov(t[p:(p+N[i]-1)], alpha, length_scale) + diag_matrix(rep_vector(sigma, N[i])));
    //print("L_K");
    //print(L_K);
    f = L_K * eta[p:(p+N[i]-1)];
    //print("f");
    //print(f);
    //print("t_e");
    //print(t_e[p:(p+N[i]-1),]);
    X = return_X(t_e[p:(p+N[i]-1),], tau);
    //print("X");
    //print(X);
    d_p = to_vector(d[p:(p+N[i] -1)]);
    y_sim[p:(p+N[i] -1)] = neg_binomial_2_log_rng(beta_0 + f + X*beta + d_p, phi);
    p = p + N[i];
  }

}

