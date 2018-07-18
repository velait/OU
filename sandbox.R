library(rstan); library(tidyverse)
seed <- 1
# all cores
options(mc.cores = parallel::detectCores())

source("OU.functions.R")

#### models ####
original_model <- stan_model("original_noncentered.stan")
original_hierarchical_model_OU <- stan_model("original_hierarchical_noncentered_OU.stan")

# The one in StanCon submission
hierarchical_noncentered <- stan_model("hierarchical_noncentered.stan")

# Gompertz removed
original_model_OU <- stan_model("original_noncentered_OU.stan")


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
long_series <- generate_a_series(intervals=1:25, mu=mu, lambda=lambda, kappa=kappa, seed=1)

long_series2 <- generate_a_series(intervals=1:25, mu=mu, lambda=lambda, kappa=kappa, seed=2)

long_series3 <- generate_a_series(intervals=1:25, mu=mu, lambda=lambda, kappa=kappa, seed=2)

long_series4 <- generate_a_series(intervals=1:25, mu=mu, lambda=lambda, kappa=kappa, seed=2)

# two long series
concatenate_two_series <- concatenate_series(list(long_series, long_series2, long_series3, long_series4))



#### Stan ####

# low number of chains for testing
chains <- 2
iter <- 4000

#### one long series, estimates OK
#### Warnings: plenty
original_long_sample <- sampling(original_model_OU,long_series, chains = chains, iter=iter)


#### two long series, same parameters

original_two_series <- sampling(original_hierarchical_model_OU, concatenate_two_series, chains = chains, iter = iter)


#### give mu and kappa to the model
# Warnings:  max TD, low BFMI, but no divergent transitions, some large Rhats.

# add kappa and mu values
concatenate_two_series_fix_par <- concatenate_two_series
concatenate_two_series_fix_par[["kappa_log"]] <- rep(log(0.1), 2)
concatenate_two_series_fix_par[["mu"]] <- rep(5, 2)

fixed_par_original_long_sample <- sampling(fixed_par_model, concatenate_two_series_fix_par, chains = chains, iter=iter)


#### TODO ####

#### simulate data with one series, show that lambda is recovered better with more data ---

# data
single_series_set <- lapply(seq(from=5, to=100, by = 5), function(x) {
  s <- generate_a_series(intervals = 1:x, mu=mu, lambda=lambda, kappa=kappa)
  s[["kappa_log"]] <- as.array(log(0.1))
  s[["mu"]] <- as.array(5)
  return(s)
})
names(single_series_set) <- as.character(seq(from=5, to=100, by = 5))

# stan samples
# single_series_set_samples <- lapply(single_series_set, function(x) sampling(fixed_par_model, x, chains=2))
# names(single_series_set_samples) <- as.character(seq(from=5, to=100, by = 5))
# 
# # save 
# save(single_series_set_samples, file="single_series_set_samples")

load(file="single_series_set_samples")

# make data frame for results
single_series_results <- matrix(NA, length(single_series_set_samples), 4)
colnames(single_series_results) <- c("length", "lower25", "mode", "upper75")

j <- 1
for(i in single_series_set_samples) {
  res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
  
  single_series_results[j, c("lower25", "mode", "upper75")] <- res
  j <- j + 1
}

single_series_results[, "length"] <- as.numeric(names(single_series_set_samples))
single_series_results <- single_series_results %>% as.data.frame()

                            
# plot
p_single_series <- plot_posteriors(single_series_results, sim_value = .5, par="Lambda ")


#### two series ----
## Data: same lambdas ----

two_series_set <- lapply(seq(from=5, to=100, by = 5), function(x) {
  
  s1 <- generate_a_series(intervals = 1:x, mu=mu, lambda=.5, kappa=kappa, seed=1)
  s1[["kappa_log"]] <- as.array(log(0.1))
  s1[["mu"]] <- as.array(5)
  
  s2 <- generate_a_series(intervals = 1:x, mu=mu, lambda=.5, kappa=kappa, seed=2)
  s2[["kappa_log"]] <- as.array(log(0.1))
  s2[["mu"]] <- as.array(5)
  
  return(concatenate_series(list(s1, s2)))
})

names(two_series_set) <- as.character(seq(from=5, to=100, by = 5))




# stan samples
# two_series_set_samples <- lapply(two_series_set, function(x) sampling(fixed_par_model, x, chains=2))
# names(two_series_set_samples) <- as.character(seq(from=5, to=100, by = 5))
# 
# save(two_series_set_samples, file="two_series_set_samples")
load(file="two_series_set_samples")

# make data frame for results
two_series_results <- matrix(NA, 2*length(two_series_set_samples), 5)
colnames(two_series_results) <- c("group", "length", "lower25", "mode", "upper75")

j <- 1
for(i in two_series_set_samples) {
  res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
  
  two_series_results[c(j, j+1), c("lower25", "mode", "upper75")] <- res
  j <- j + 2
}


two_series_results[, "length"] <- rep(as.numeric(names(two_series_set_samples)), each=2)
two_series_results <- two_series_results %>% as.data.frame()
two_series_results[, "group"] <- rep(c("A","B"), length(two_series_set_samples))


# plot
p_two_series_plot <- ggplot(two_series_results, aes(x=length, y=mode, group=group, color=group)) + 
  geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
  geom_line() +
  geom_point()  + geom_hline(yintercept = c(.5), linetype="dashed", color="black") + labs(y="Estimate", x="Length", title="Lambda estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, 7))









## Data: different lambdas: .3 and .7 ----
two_diff_series_set <- lapply(seq(from=5, to=100, by = 5), function(x) {
  
  s1 <- generate_a_series(intervals = 1:x, mu=mu, lambda=.3, kappa=kappa, seed=1)
  s1[["kappa_log"]] <- as.array(log(0.1))
  s1[["mu"]] <- as.array(5)
  
  s2 <- generate_a_series(intervals = 1:x, mu=mu, lambda=.7, kappa=kappa, seed=2)
  s2[["kappa_log"]] <- as.array(log(0.1))
  s2[["mu"]] <- as.array(5)
  
  return(concatenate_series(list(s1, s2)))
})

names(two_diff_series_set) <- as.character(seq(from=5, to=100, by = 5))




# stan samples
# two_diff_series_set_samples <- lapply(two_diff_series_set, function(x) sampling(fixed_par_model, x, chains=2))
# names(two_diff_series_set_samples) <- as.character(seq(from=5, to=100, by = 5))
# 
# # save
# save(two_diff_series_set_samples, file="two_diff_series_set_samples")
load(file="two_diff_series_set_samples")

# make data frame for results
two_diff_series_results <- matrix(NA, 2*length(two_diff_series_set_samples), 5)
colnames(two_diff_series_results) <- c("group", "length", "lower25", "mode", "upper75")

j <- 1
for(i in two_diff_series_set_samples) {
  res <- summary(i)$summary[grep("lambda\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
  
  two_diff_series_results[c(j, j+1), c("lower25", "mode", "upper75")] <- res
  j <- j + 2
}


two_diff_series_results[, "length"] <- rep(as.numeric(names(two_diff_series_set_samples)), each=2)
two_diff_series_results <- two_diff_series_results %>% as.data.frame()
two_diff_series_results[, "group"] <- rep(c("A","B"), length(two_diff_series_set_samples))


# plot
two_diff_series_plot <- ggplot(two_diff_series_results, aes(x=length, y=mode, group=group, color=group)) + 
  geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
  geom_line() +
  geom_point()  +
  geom_hline(yintercept = c(.3), linetype="dashed", color="#f1a340") + 
  geom_hline(yintercept = c(.7), linetype="dashed", color="#998ec3") + 
  labs(y="Estimate", x="Length", title="Lambda estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, 10))




## Same hierarchical models but separately for the data sets to assess the improvement over non hierarchical model

## lambda = .3 ----
single_series_set2 <- lapply(seq(from=5, to=100, by = 5), function(x) {
  s <- generate_a_series(intervals = 1:x, mu=mu, lambda=.3, kappa=kappa)
  s[["kappa_log"]] <- as.array(log(0.1))
  s[["mu"]] <- as.array(5)
  return(s)
})
names(single_series_set2) <- as.character(seq(from=5, to=100, by = 5))

# # stan samples
# single_series_set_samples2 <- lapply(single_series_set2, function(x) sampling(fixed_par_model, x, chains=2))
# names(single_series_set_samples2) <- as.character(seq(from=5, to=100, by = 5))
# 
# # save
# save(single_series_set_samples2, file="single_series_set_samples2")
load(file="single_series_set_samples2")

# get df with mode and 50% intervals
single_series2_results <- samples_df(single_series_set_samples2)

# plot posteriors w simulation value
p_single_series2 <- plot_posteriors(single_series2_results, sim_value = .3, par="Lambda ")



## lambda = .7 ----

single_series_set3 <- lapply(seq(from=5, to=100, by = 5), function(x) {
  s <- generate_a_series(intervals = 1:x, mu=mu, lambda=.7, kappa=kappa)
  s[["kappa_log"]] <- as.array(log(0.1))
  s[["mu"]] <- as.array(5)
  return(s)
})
names(single_series_set3) <- as.character(seq(from=5, to=100, by = 5))

# stan samples
# single_series_set_samples3 <- lapply(single_series_set3, function(x) sampling(fixed_par_model, x, chains=2))
# names(single_series_set_samples3) <- as.character(seq(from=5, to=100, by = 5))

# save
# save(single_series_set_samples3, file="single_series_set_samples3")

# load
load(file="single_series_set_samples3")

# get df with mode and 50% intervals
single_series3_results <- samples_df(single_series_set_samples3)

# plot posteriors w simulation value
p_single_series3 <- plot_posteriors(single_series3_results, sim_value = .7, par="lambda")

## difference plots ----

# Does the .3 graph from the hierarchical model converge faster compared to the single series model?

df_0.3 <- rbind(cbind(single_series2_results, group="single"), two_diff_series_results %>% filter(group=="A") %>% mutate(group="hierarchical"))

p_diff_0.3 <- plot_hier_posteriors(df_0.3, sim_value = .3) + labs(title="lambda = 0.3")


# lambda = .7

df_0.7 <- rbind(cbind(single_series3_results, group="single"), two_diff_series_results %>% filter(group=="B") %>% mutate(group="hierarchical"))

p_diff_0.7 <- plot_hier_posteriors(df_0.7, sim_value = .7) + labs(title="lambda = 0.7")


# lambda = .5
# df_0.5 <- rbind(two_series_results, )

p_diff_0.5 <- ggplot(rbind(single_series_results %>% mutate(group="single"), two_series_results), aes(x=length, y=mode, group=group, color=group)) + 
  geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
  geom_line() +
  geom_point()  + geom_hline(yintercept = .5, linetype="dashed", color="black") + labs(y="Estimate", x="Length", title="lambda = 0.5") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3", "black")) + scale_y_continuous(limits = c(0, 7))



#### find representative plots ####



x1 <- generate_a_series(intervals = 1:100, mu=mu, lambda=.1, kappa=5)
x2 <- generate_a_series(intervals = 1:100, mu=mu, lambda=.5, kappa=1)
x3 <- generate_a_series(intervals = 1:100, mu=mu, lambda=.9, kappa=0.5556)

par(mfrow=c(3,1)) 
plot(y= x1[["observations"]], x=1:100, type="l", ylim=c(0, 10), main="lambda=0.1, constant sigma")
plot(y= x2[["observations"]], x=1:100, type="l", ylim=c(0, 10), main="lambda=0.5, constant sigma")
plot(y= x3[["observations"]], x=1:100, type="l", ylim=c(0, 10), main="lambda=0.9, constant sigma")




#### SIGMA ----

### SINGLE SERIES

sigma_series <- lapply(seq(from=5, to=100, by = 5), function(x) {
  s <- generate_a_sigma_series(intervals = 1:x, mu=mu, lambda=.5, kappa=kappa)
  s[["kappa_log"]] <- as.array(log(0.1))
  s[["mu"]] <- as.array(5)
  return(s)
})



# test <- lapply(sigma_series, function(x) sampling(original_model_OU, x, chains=2, iter=1000))
# 
# save(test, file="sigma_series_test")
load(file="sigma_series_test")

names(sigma_series) <- as.character(seq(from=5, to=100, by = 5))
names(test) <- names(sigma_series)


test_plot_sigma <- plot_posteriors(samples_df_par(test, par="sigma"), sim_value = .5, par="Sigma ")
test_plot_kappa <- plot_posteriors(samples_df_par(test, par='kappa'), sim_value = 0.5, par="Kappa ")





#### HIERARCHICAL, SIGMA = 0.5

sigma_set <- lapply(seq(from=5, to=100, by = 5), function(x) {
  
  s1 <- generate_a_sigma_series(intervals = 1:x, mu=mu, lambda=.5, kappa=kappa, seed=1)
  s1[["kappa_log"]] <- as.array(log(0.1))
  s1[["mu"]] <- as.array(5)
  
  s2 <- generate_a_sigma_series(intervals = 1:x, mu=mu, lambda=.5, kappa=kappa, seed=2)
  s2[["kappa_log"]] <- as.array(log(0.1))
  s2[["mu"]] <- as.array(5)
  
  return(concatenate_series(list(s1, s2)))
})


# hier_sigma_test <- lapply(sigma_set, function(x) sampling(original_hierarchical_model_OU, x, chains=2, iter=1000))

names(hier_sigma_test) <- as.character(seq(from=5, to=100, by = 5))

# save(hier_sigma_test, file="hier_sigma_test")

load(file="hier_sigma_test")



# make data frame for results
two_sigma <- matrix(NA, 2*length(hier_sigma_test), 5)
colnames(two_sigma) <- c("group", "length", "lower25", "mode", "upper75")

j <- 1
for(i in hier_sigma_test) {
  res <- summary(i)$summary[grep("sigma\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
  
  two_sigma[c(j, j+1), c("lower25", "mode", "upper75")] <- res
  j <- j + 2
}


two_sigma[, "length"] <- rep(as.numeric(names(hier_sigma_test)), each=2)
two_sigma <- two_sigma %>% as.data.frame()
two_sigma[, "group"] <- rep(c("A","B"), length(hier_sigma_test))


# plot
two_sigma_plot <- ggplot(two_sigma, aes(x=length, y=mode, group=group, color=group)) + 
  geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
  geom_line() +
  geom_point()  +
  geom_hline(yintercept = c(.5), linetype="dashed", color="black") + 
  labs(y="Estimate", x="Length", title="Sigma estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, 7))






#### KAPPA
# make data frame for results
two_kappa <- matrix(NA, 2*length(hier_sigma_test), 5)
colnames(two_kappa) <- c("group", "length", "lower25", "mode", "upper75")

j <- 1
for(i in hier_sigma_test) {
  res <- summary(i)$summary[grep("kappa\\[", rownames(summary(i)$summary)),c("25%", "50%", "75%")]
  
  two_kappa[c(j, j+1), c("lower25", "mode", "upper75")] <- res
  j <- j + 2
}


two_kappa[, "length"] <- rep(as.numeric(names(hier_sigma_test)), each=2)
two_kappa <- two_kappa %>% as.data.frame()
two_kappa[, "group"] <- rep(c("A","B"), length(hier_sigma_test))


# plot
two_kappa_plot <- ggplot(two_kappa, aes(x=length, y=mode, group=group, color=group)) + 
  geom_errorbar(aes(ymin=lower25, ymax=upper75), width=1) +
  geom_line() +
  geom_point()  +
  geom_hline(yintercept = c(.5), linetype="dashed", color="black") + 
  labs(y="Estimate", x="Length", title="Kappa estimates with 50% error bars vs. series length") + theme_bw() + scale_color_manual(values=c("#f1a340", "#998ec3")) + scale_y_continuous(limits = c(0, 5))



#### CHECK SIGMA WITH MULTIPLE SERIES

sigma_list <- concatenate_series(lapply(1:5, function(i) generate_a_sigma_series(intervals = 1:25, mu=mu, lambda=.5, kappa=kappa, seed=i)))

sigma_hier_samples <- sampling(original_hierarchical_model_OU, sigma_list, chains=chains, iter=2000)

