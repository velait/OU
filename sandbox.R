seed <- 1
source("OU.functions.R")

#### models ####
original_model <- stan_model("original_noncentered.stan")
original_hierarchical_model <- stan_model("original_hierarchical_noncentered.stan")
hierarchical_noncentered <- stan_model("hierarchical_noncentered.stan")


#### data ####

# original model

kappa <- .1
lambda <- .5
mu <- 8
t.df <- 5

long_series <- generate_data_original(kappa = kappa, lambda = lambda, mu = mu, t.df = t.df, intervals = 1:25, seed = seed)

