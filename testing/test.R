library(tibble)
library(dplyr)
library(tidyverse)
library(magrittr)
library(plotly)
theme_set(theme_bw(20))
#### Pooled vs. partially pooled vs. non pooled models

compare_n_series <- 20
chains <- 3
iter <- 2e3
intervals <- 1:30

mu_val <- rep(5, compare_n_series)
sigma_val <- rep(0.2, compare_n_series)
lambda_val <- oup_invG_lambda(compare_n_series, 4, 0.4)

diff_compare_data <- generate_student_set(n_series = compare_n_series,
                                          student_df = 7,
					  mu = mu_val,
					  sigma =  sigma_val,
					  lambda = lambda_val,
					  intervals = intervals)

#### MODELS ####
hierarchical_student_t_oup <- stan_model("hierarchical_student_t_oup.stan")

#### SAMPLES ####
diff_partially_pooled_samples <- sampling(hierarchical_student_t_oup, diff_compare_data, iter=iter, chains=chains, init=0.5)

#### RESULTS ####

# Get hyper parameter summaries in list
diff_compare_samples_list <- list(partially_pooled_samples=diff_partially_pooled_samples) 

diff_compare_summary <- list()
for(x in names(diff_compare_samples_list)) {
  diff_compare_summary[[x]] <- summary(diff_compare_samples_list[[x]])$summary[grep("lambda\\[|sigma\\[|mu\\[", rownames(summary(diff_compare_samples_list[[x]])$summary)), c("25%", "50%", "75%")] %>%
    as.data.frame() %>%
    rownames_to_column(var="parameter") %>%
    set_colnames(c("parameter", "lower", "mode", "upper"))
}

#### PLOT ALL ####

par <- "lambda"
ord <- diff_compare_data[[paste0(par, "_values")]] %>% order
partial_estimate <- diff_compare_summary[["partially_pooled_samples"]] %>% 
                         filter(str_detect(parameter, par)) %>% 
			 mutate(series = str_extract(parameter, "\\d+") %>% 
			 as.numeric(), 
			 model="Partially pooled")

m <- max(c(lambda_val, partial_estimate$mode))
plot(lambda_val, partial_estimate$mode, main = paste(par, "; r=", round(cor(lambda_val, partial_estimate$mode), 2)), xlim = c(0, m), ylim = c(0,m))

abline(0,1, lty = 2)






