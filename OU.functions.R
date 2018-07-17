# Functions

# Alternative generator
ou_simulator <- function (T, mu, lambda, kappa, x0 = NULL, seed=1) {
  seed <- seed
  x <- c()
  
  # Initial value
  if (is.null(x0)) {
    x[[1]] <- mu + rnorm(1) * sqrt(kappa)
  } else {
    x[[1]] <- x0
  }
  
  # Consecutive values
  for (t in 2:T) {
    x[[t]] <- mu - (mu - x[[t-1]]) * exp(-lambda) + rnorm(1) * sqrt(kappa * (1 - exp(-2 * lambda)))
  }
  
  #gompertz assumptions, poisson model
  list(observations = rpois(T,exp(x)), T=T, n_series=1, samples_per_series=as.array(T), time=1:T)
  
  
}




# original generator
generateStanData <- function(kappa,
                             lambda,
                             mu,
                             intervals,
                             t.df = Inf,
                             seed = 1){
  set.seed(seed)
  N <- length(intervals)
  lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  out.data$latent_value <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  
  out.data$value <- rpois(N,exp(out.data$latent_value))
  
  
  out.data$time <- intervals
  out.data$replicates <- 1L
  out.data$replicate_samples <- array(length(intervals))
  out.data$NSMPL <- length(intervals)
  out.data
}

# original generator w different parameter names
generate_a_series <- function(kappa,
         lambda,
         mu,
         intervals,
         t.df = Inf,
         seed = 1){
  set.seed(seed)
  N <- length(intervals)
  lv.variates <- rnorm(N)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=N) else 1
  out.data <- list()
  # out.data$latent_value <- as.vector(t(L) %*% (rnorm(N) * scale)) + mu
  out.data$observations <- rpois(N,exp(as.vector(t(L) %*% (rnorm(N) * scale))+mu))
  out.data$time <- intervals
  out.data$n_series <- 1L
  out.data$samples_per_series <- array(length(intervals))
  out.data$T <- length(intervals)
  out.data
}

# Generate one long for the original
generate_data_original <- function(kappa=0.1, lambda=0.5, mu=5, intervals, t.df = Inf, seed = 1){
  set.seed(seed)
  n_obs <- length(intervals)
  lv.variates <- rnorm(n_obs)
  dt <- outer(intervals,intervals,function(x,y) abs(x-y))
  x <- kappa * exp(-lambda*dt)
  L <- chol(x)
  
  scale <- if(is.finite(t.df)) rep(sqrt(rgamma(1,t.df/2,(t.df-2)/2)),each=n_obs) else 1
  
  data <- list()
  Y <- matrix(NA, nrow=n_series, ncol=length(intervals))
  
  
    X_latent <- as.vector(t(L) %*% (rnorm(n_obs) * scale)) + mu
    Y <- rpois(n_obs,exp(X_latent))
  
  
  data[["n_series"]] <- 1           # number of time series
  data[["samples_per_series"]] <- c(length(intervals)) # number of time series
  data[["T"]] <- length(intervals)  # number of time points
  data[["observations"]] <- Y                  # observations
  data[["time"]] <- intervals       # observation times
  
  return(data)
}

# Concatenate several short series

generate_data_set2 <- function(n_series, kappa, lambda, mu, intervals, t.df = Inf, seed = 1){
  
  setlist <- list()
  
  for(i in 1:n_series) {
    setlist[[i]] <- generateStanData(kappa,
                                     lambda,
                                     mu,
                                     intervals,
                                     t.df = Inf,
                                     seed = i)
  }
  
  T <- n_series*length(intervals)
  samples_per_series <- rep(length(intervals), n_series)
  observations <- c()
  time <- c()
  for(i in 1:length(setlist)) {
    observations <- c(observations, setlist[[i]][["value"]])
    time <- c(time, setlist[[i]][["time"]])
  }
  
  
  return(list(T=T, samples_per_series=samples_per_series, observations=observations, time=time, n_series=n_series))
  
}

# Generate a data set
generate_data_set <- function(kappa=0.1, lambda=0.5, mu=5, intervals, n_series,
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


# concatenate series
concatenate_series <- function(s1, s2) {
  
  obs <- c(s1[["observations"]], s2[["observations"]])
  t <- c(s1[["time"]], s2[["time"]])
  n <- 2
  s <- c(s1[["samples_per_series"]], s2[["samples_per_series"]])
  T <- length(obs)
  
  
  list(observations=obs, time=t, n_series=n, samples_per_series=s, T=T)
  
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




