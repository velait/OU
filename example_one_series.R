seed <- 1
library(rstan)
source("OU.functions.R")

# all cores
options(mc.cores = parallel::detectCores())

# model 
original_model <- stan_model("original_noncentered.stan")


#### data ####
#model parameters
kappa <- .1
lambda <- 1
mu <- 10
t.df <- 5
intervals <- 1:25

long_series <- generate_a_series(kappa, lambda, mu, intervals = intervals, t.df = t.df, seed = seed)

#### Stan ####

# low number of chains for testing
chains <- 2
iter <- 2000

# one long series, estimates OK
# Warnings: BFMI low, max treedepth exceedings
original_long_sample <- sampling(original_model,long_series, chains = chains, iter=iter)