# Functions

# Generate a data set
generate_data_set <- function(kappa=0.1, lambda=1, mu=5, intervals, n_series,
                              t.df = Inf, seed = 1){
  
  set.seed(seed)
  n_obs <- length(intervals)
  lv.variates <- rnorm(n_obs)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=n_obs) else 1
  
  data <- list()
  Y <- matrix(NA, nrow=n_series, ncol=length(intervals))
  
  for(i in 1:n_series) {
    X_latent <- as.vector(t(L) %*% (rnorm(n_obs) * scale)) + mu
    Y[i, ] <- rpois(n_obs,exp(X_latent))
  }
  
  data[["N"]] <- n_series           # number of time series
  data[["T"]] <- length(intervals)  # number of time points
  data[["Y"]] <- Y                  # observations
  data[["time"]] <- intervals       # observation times
  
  return(data)
}

# HPDI
success_rate_HPDI <- function(stan_fit, parameter="lambda", real_value=1, alpha=0.95) {
  
  pos <- rstan::extract(stan_fit)[[parameter]]
  
  success <- 0
  for(i in 1:ncol(pos)) {
    hpdi <- HPDI(pos[,i], prob=alpha)
    
    if(hpdi[1] < real_value && real_value < hpdi[2]) {
      success <- success + 1
    }
    
  }
  success/ncol(pos)
}

# Percentiles
success_rate_quantiles <- function(stan_fit, parameter="lambda", real_value=1, success_limit=0.95) {
  
  pos <- rstan::extract(stan_fit)[[parameter]]
  
  success <- 0
  for(i in 1:ncol(pos)) {
    q <- quantile(pos[,i], probs = c(1-success_limit, 0.5, success_limit))
    
    if(q[1] < real_value && real_value < q[3]) {
      success <- success + 1
    }
    
  }
  success/ncol(pos)
}
