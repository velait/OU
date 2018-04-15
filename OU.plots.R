# Plot results
success_plots <- list()
for(p in c("kappa", "lambda", "mu")) {
  
  for(t in c("quan_table", "hpdi_table")) {
    
    tbl <- as.data.frame(results[[p]][[t]])
    tbl$n_series <- as.factor(row.names(tbl))
    tbl1 <- tbl[1:5,]
    tbl1_melt <- melt(tbl1)
    
    success_plots[[p]][[t]][[1]] <- ggplot(data=tbl1_melt, aes(x=variable, y=value, colour=n_series, group=n_series)) + geom_line() + coord_cartesian(ylim=c(0,1)) + labs(x="observations", y="proportion", subtitle=paste(p,ifelse(t=="quan_table", "(percentiles)", "(HPDI)")))
    
    # 5 obs, series 5 -> 100
    tbl2 <- data.frame(prop=tbl[,"5"], n_series=as.numeric(levels(tbl$n_series))[tbl$n_series])
    
    success_plots[[p]][[t]][[2]] <- ggplot(data=tbl2, aes(x=n_series, y=prop)) + geom_line() + labs(x="series", y="proportion", subtitle=paste(p,ifelse(t=="quan_table", "(percentiles)", "(HPDI)"))) + coord_cartesian(ylim=c(0,1))
    
  }
}


# running times
running_time <- matrix(NA, nrow=length(n_series), ncol=length(n_obs))
colnames(running_time) <- as.character(n_obs)
rownames(running_time) <- as.character(n_series)
for(i in 1:length(fit_list)) {
  time <- sapply(fit_list[[i]]@sim[[1]],function(x) attr(x, "elapsed_time" )) %>% sum
  
  obs <- sub(".*\\b(\\d+).*", "\\1", names(quan)[[i]])
  series <- sub("\\D*(\\d+).*", "\\1", names(quan)[[i]])
  
  running_time[series, obs] <- time
}

running_time <- running_time %>% round()


#sampling times in minutes for 5 observations
running_time_long <- data.frame(minutes=round(running_time[,"5"]/(60)),
                                n_series=rownames(running_time)) %>%  kable()
