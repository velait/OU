diff_compare_plots <- list()
diff_partially_pooled_errorbar_plots <- list()
diff_no_pooling_errorbar_plots <- list()

for(par in c("lambda", "sigma", "mu")) {
  
  #ord <- diff_compare_data[[paste0(par, "_values")]] %>% order
  ord <- 1:length(diff_compare_data[[paste0(par, "_values")]])

  # simulation_value <- oup_simulation_parameters[par]
  
  partial_estimate <- diff_compare_summary[["partially_pooled_samples"]] %>% 
                         filter(str_detect(parameter, par)) %>% 
			 mutate(series = str_extract(parameter, "\\d+") %>% 
			 as.numeric(), 
			 model="Partially pooled")
  
  partial_estimate <- partial_estimate[match(ord, partial_estimate$series),] %>% 
                                          mutate(ord = 1:compare_n_series)
  
  estimates <- partial_estimate %>% dplyr::select(-parameter)
  
  p <- ggplot() +
    geom_point(data=estimates, aes(x=ord, y=mode, shape=model))  +
    labs(x="Series", y="Posterior estimate", title=paste0("\\", par)) + 
    guides(shape=guide_legend("")) +
    # geom_hline(yintercept=pooled_estimate, linetype = "dashed") +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]] %>%
    sort() %>%
    as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16))
  
  diff_compare_plots[[par]] <- p
  
  q <- ggplot()  + 
    labs(x="Series", y="Posterior estimate", title=par) +
    guides(shape=guide_legend("")) +
    geom_line(data=diff_compare_data[[paste0(par, "_values")]]%>%
    sort() %>%
    as_tibble(), aes(y=value, x=1:compare_n_series)) +
    scale_shape_manual(values=c(1, 16)) + 
    geom_errorbar(data=estimates %>% filter(model=="Partially pooled"), aes(x=ord, ymin=lower, ymax=upper)) 
  
  diff_partially_pooled_errorbar_plots[[par]] <- q
  
}


