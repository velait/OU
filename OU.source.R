# Source for the objects needed in the Rmd. 

#### Parameters ####
sigma <- 0.1
lambda <- 0.1
mu <- 5
t.df <- 7

#### Example plots ####


# lambda = 0.01, 0.1, 1

oup_example_plot <- sapply(c(0.01, 0.1, 1), function(x) generate_a_series(sigma = sigma, lambda = x, mu = mu, intervals = 1:200, t.df=Inf, seed = 1)[["observations"]]) %>% as_tibble() %>% set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value, color=variable)) + geom_line() + scale_y_continuous(limits = c(4, 7)) +  theme_bw() + labs(x="Time", y="Value", title="mu = 5, sigma = 0.1") + scale_color_manual(values=c('#4daf4a','#377eb8', '#e41a1c')) + theme(legend.title=element_blank())


# sigma = 0.01, 0.1, 1
oup_example_plot2 <- sapply(c(0.05, 0.1, 0.2), function(x) generate_a_series(sigma = x, lambda = lambda, mu = mu, intervals = 1:200, t.df=Inf, seed = 1)[["observations"]]) %>% as_tibble() %>% set_colnames(c("sigma = 0.05", "sigma = 0.1", "sigma = 0.2")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value, color=variable)) + geom_line() + theme_bw() + labs(x="Time", y="Value", title="mu = 5, lambda = 0.1") +  scale_color_manual(values=c('#4daf4a','#377eb8', '#e41a1c')) + theme(legend.title=element_blank())

####  ####