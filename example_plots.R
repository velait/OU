# Generate representative examples of series 
# Set sigma low to see mean-reversion

library(reshape2)
library(magrittr)
sigma <- 0.1
mu <- 5
t.df <- 7

# lambda = 0.01, 0.1, 1

oup_example_plot <- sapply(c(0.01, 0.1, 1), function(x) generate_a_series(sigma = 0.1, lambda = x, mu = mu, intervals = 1:200, t.df=Inf, seed = 1)[["observations"]]) %>% as_tibble() %>% set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value)) + geom_line() + scale_y_continuous(limits = c(2.5, 7.5)) + facet_grid(cols = vars(variable)) + theme_bw() + theme(legend.position="none") + labs(x="Time", y="Value", title="mu = 5, sigma = 0.1")

jpeg("oup_example_plot.jpg", width = 750, height = 300)
oup_example_plot
dev.off()


oup_example_plot2 <- sapply(c(0.01, 0.1, 1), function(x) generate_a_series(sigma = 0.05, lambda = x, mu = mu, intervals = 1:200, t.df=Inf, seed = 1)[["observations"]]) %>% as_tibble() %>% set_colnames(c("lambda = 0.01", "lambda = 0.1", "lambda = 1")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value)) + geom_line() + scale_y_continuous(limits = c(2.5, 7.5)) + facet_grid(cols = vars(variable)) + theme_bw() + theme(legend.position="none") + labs(x="Time", y="Value", title="mu = 5, sigma = 0.05")

jpeg("oup_example_plot2.jpg", width = 750, height = 300)
oup_example_plot2
dev.off()



student_oup_example_plot <- sapply(c(3, 7, Inf), function(x) generate_a_series(sigma = 0.05, lambda = 0.1, mu = mu, intervals = 1:200, t.df=x, seed = 1)[["observations"]]) %>% as_tibble() %>% set_colnames(c("df = 3", "df = 10", "df = Inf")) %>%  mutate(time=1:200) %>% melt(id.vars="time") %>% ggplot(aes(x=time, y=value)) + geom_line() + scale_y_continuous(limits = c(3, 7)) + facet_grid(cols = vars(variable)) + theme_bw() + theme(legend.position="none") + labs(x="Time", y="Value", title="mu = 5, sigma = 0.05")

jpeg("student_oup_example_plot.jpg", width = 750, height = 300)
student_oup_example_plot
dev.off()