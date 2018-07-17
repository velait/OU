//Copyright 2017 Aaron Goodman <aaronjg@stanford.edu>. Licensed under the GPLv3 or later.
data{
  int <lower=0> T;               // number of samples
  int <lower=0> n_series;        // number of series
  int <lower=0> samples_per_series[n_series]; 
  int <lower=0> observations[T]; // observation concatenated
  vector [T] time;               // observation times concatenated
}
transformed data{
  
  vector[T] observation_vec = to_vector(observations);
  vector[T] time_vec = to_vector(time);
  
  int sample_start[n_series];
  sample_start[1] = 1;
  if(n_series > 1) {
    for(i in 2:n_series)
      sample_start[i] = sample_start[i-1] + samples_per_series[i-1];
  }
  
}

parameters{
  real lambda_log;
  real kappa_log;
  real mu;
  real <lower=2> student_df;
  vector[T] epsilon; // t-process observations
}

transformed parameters{
  real sigma_log = 0.5*(kappa_log + lambda_log + log2());
  real sigma = exp(sigma_log);
  real lambda = exp(lambda_log);
  real kappa = exp(kappa_log);
  real kappa_inv = exp(-kappa_log);
  real kappa_sqrt = exp(0.5*kappa_log);
  vector[T] latent_observation;
  
  // For each raw latent observation, transform it using the
  // conditional expression for the OU process.
  for(i in 1:n_series){
    
    int n = samples_per_series[i];
    int offset = sample_start[i];
    vector [n-1] time_diff = segment(time_vec, offset + 1, n - 1) -
      segment(time_vec,offset,n-1);
    real cum_squares = 0;
    real last_t = -1;
    real lv;
    for(k in 0:n-1){
      real lv_raw = epsilon[offset + k];
      if(k == 0){
        //For the first latent observation use the stationary distribution.
        lv = mu + lv_raw * kappa_sqrt;
      }else{
        real t = time_diff[k];
        real exp_neg_lambda_t = exp(-t*lambda);
        real sd_scale = kappa_sqrt .* sqrt(1-square(exp_neg_lambda_t));
        lv = mu - (mu - lv) .* exp_neg_lambda_t + lv_raw .* sd_scale;
        last_t = t;
      }
      latent_observation[offset+k] = lv;
    }
    
  }
}
model {
  target += lambda_log;
  target += kappa_log;
  student_df ~ gamma(2,.1);
  
  // Increment the log probability according to the conditional expression for
  // the unit multivariate t-distribution.
  for(i in 1:n_series){
    int n = samples_per_series[i];
    int offset = sample_start[i];
    vector [n] sq_lv = square(segment(epsilon,offset,n));
    vector [n] cum_squares = cumulative_sum(append_row(0,sq_lv[1:n-1]));    
    
    for(k in 1:T){
      target +=  (lgamma((student_df + k) * 0.5) - lgamma((student_df+ k - 1 )* 0.5));     
      target += -0.5 * (student_df + k) * log1p(sq_lv[k] / (student_df + cum_squares[k] - 2));
      target += -0.5 * log(student_df + cum_squares[k] - 2);
    }
  }
  
  // Add the log probability for the observations given the latent observations
  observations ~ poisson_log(latent_observation);
  
  // Prior probabilities (not the original ones).
  lambda ~ normal(0,5);
  kappa  ~ normal(0,5);
  mu     ~ normal(0,5);
}
