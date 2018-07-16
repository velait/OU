seed <- 1
source("OU.functions.R")

#### models ####
original_model <- stan_model("original_noncentered.stan")
original_hierarchical_model <- stan_model("original_hierarchical_noncentered.stan")

# The one in StanCon submission
hierarchical_noncentered <- stan_model("hierarchical_noncentered.stan")


#### data ####
#model parameters
kappa <- .1
lambda <- .5
mu <- 8
t.df <- 5

# original model ----
test_data_original <- generateStanData(kappa, lambda, mu, intervals = 1:50, t.df = Inf, seed = seed)



# data with different variable names
long_series <- generate_one_series(kappa, lambda, mu, intervals = 1:50, t.df = Inf, seed = seed)

long_series2 <- generate_one_series(kappa, lambda, mu, intervals = 1:50, t.df = Inf, seed = 2)

# two long series
concatenate_two_series <- concatenate_series(long_series, long_series2)



#### Stan ####

# low number of chains
chains <- 2
iter <- 4000

# one long series, estimastes OK
# Warnings: divergent transitions, BFMI low, max treedepth exceedings
original_long_sample <- sampling(original_model,long_series, chains = chains, iter=iter)


# two long series, same parameters
original_two_series <- sampling(original_hierarchical_model, concatenate_two_series, chains = chains, iter = iter)
