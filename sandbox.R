seed <- 1
source("OU.functions.R")

#### models ####
original_model <- stan_model("original_noncentered.stan")
original_hierarchical_model <- stan_model("original_hierarchical_noncentered.stan")

# The one in StanCon submission
hierarchical_noncentered <- stan_model("hierarchical_noncentered.stan")

# lambda and kappa fixed
fixed_par_model <- stan_model("fixed_par_original_hierarchical_noncentered.stan")

#### data ####
#model parameters
kappa <- 0.1
lambda <- .5
mu <- 5
t.df <- 5
intervals <- 1:100

# data with different variable names
long_series <- ou_simulator(T=100, mu=mu, lambda=lambda, kappa=kappa, x0=6, seed=1)

long_series2 <- ou_simulator(T=100, mu=mu, lambda=lambda, kappa=kappa, x0=6, seed=2)

# two long series
concatenate_two_series <- concatenate_series(long_series, long_series2)



#### Stan ####

# low number of chains for testing
chains <- 2
iter <- 4000

#### one long series, estimates OK for lambda and mu
original_long_sample <- sampling(original_model,long_series, chains = chains, iter=iter)


#### two long series, same parameters

original_two_series <- sampling(original_hierarchical_model, concatenate_two_series, chains = chains, iter = iter)


#### give mu and kappa to the model
# Warnings:  max TD, low BFMI, but no divergent transitions, some large Rhats.

# add kappa and mu values
concatenate_two_series_fix_par <- concatenate_two_series
concatenate_two_series_fix_par[["kappa_log"]] <- rep(log(0.1), 2)
concatenate_two_series_fix_par[["mu"]] <- rep(5, 2)

fixed_par_original_long_sample <- sampling(fixed_par_model, concatenate_two_series_fix_par, chains = chains, iter=iter)


