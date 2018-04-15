// Hierarchical model with equal time points, non-centered parametization
data{
  int <lower=0> N;          //Number of time series
  int <lower=0> T;          //Number of time points, same in each
  int<lower=0> Y[N,T];      //observations matrix as 2D array
  vector[T] time;           //observation times
}
transformed data{
  vector[T] time_vec = to_vector(time);
}
parameters{
  vector[N] lambda_log; 
  vector[N] kappa_log;  
  vector[N] mu; 
  real <lower=2> student_df; 
  matrix[N,T] error_terms;   //Error terms
}
transformed parameters{
  vector[N] sigma_log = 0.5*(kappa_log + lambda_log + log2());
  vector<lower=0>[N] sigma = exp(sigma_log);
  vector<lower=0>[N] lambda = exp(lambda_log);
  vector<lower=0>[N] kappa = exp(kappa_log);
  vector<lower=0>[N] kappa_inv = exp(-kappa_log);
  vector<lower=0>[N] kappa_sqrt = exp(0.5*kappa_log);
  matrix[N, T] X_latent; 
  
  // Relate error terms to latent values
  for(i in 1:N){
    
    vector [T-1] delta_t = segment(time_vec, 2, T - 1) - segment(time_vec,1,T-1);
    real cum_squares = 0;
    real X;
    
    for(k in 1:T){
      real epsilon = error_terms[i,k];
      if(k == 1){
        //For the first latent value use the stationary distribution.
        X = mu[i] + epsilon * kappa_sqrt[i];
      }else{
        real t = delta_t[k-1];
        real exp_neg_lambda_t = exp(-t*lambda[i]);
        real sd_scale = kappa_sqrt[i] .* sqrt(1-square(exp_neg_lambda_t));
        X = mu[i] - (mu[i] - X) .* exp_neg_lambda_t + epsilon .* sd_scale;
      }
      X_latent[i, k] = X;
    }
  }
  
}
model{
  target += lambda_log;
  target += kappa_log;
  student_df ~ gamma(2,.1);
  
  // Increment the log probability according to the conditional expression for
  // the error terms 
  for(i in 1:N){
    // row vector to regular
    vector[T] error_terms2 = to_vector(square(segment(error_terms[i],1,T)));
    vector[T] cum_error_terms2 = cumulative_sum(append_row(0,error_terms2[1:T-1]));    
    
    for(k in 1:T){
      target +=  (lgamma((student_df + k) * 0.5) - lgamma((student_df+ k - 1 )* 0.5));     
      target += -0.5 * (student_df + k) * log1p(error_terms2[k] / (student_df + cum_error_terms2[k] - 2));
      target += -0.5 * log(student_df + cum_error_terms2[k] - 2);
    }
  }
  // Log likelihood for the observations given the latent values
  for(i in 1:N) {
    Y[i] ~ poisson_log(X_latent[i]);
  }
  
  // Prior probabilities.
  lambda ~ gamma(2,2);
  kappa ~ gamma(2,2);
  mu ~ normal(0,5); 
  
}
