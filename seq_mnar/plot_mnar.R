

process_df <- function(label, name="SeqMNAR_continuous"){
  
  direct_1 = "coeff_uniform_neg1_1"
  direct_2 = "coeff_uniform_neg5_15"
  direct_3 = "coeff_uniform_0_2"
  
  df_1 = read.csv(paste0(direct_1, "/", name, ".csv"), sep = ",")
  df_2 = read.csv(paste0(direct_2, "/", name, ".csv"), sep = ",")
  df_3 = read.csv(paste0(direct_3, "/", name, ".csv"), sep = ",")
  
  df = data.frame("size" = df_1$n, "low"= df_1$accept_rate, "medium"= df_2$accept_rate, "high"= df_3$accept_rate)

  plot_df = melt(df, id.vars="size", variable.name = "model", value.name = "rate") 
  if (is.null(label) == F) plot_df$data_type = label
  
  return(plot_df)
}

plot_func <- function(plot_d, title){
  
  p <- ggplot(plot_d, aes(size, rate, col=model)) + 
    geom_line() +  
    geom_point() + 
    geom_hline(yintercept=0, linetype="dashed", color = "black") + 
    xlim(500, 15500) + 
    ylim(0, 1) + 
    labs(x="Sample size", y="Acceptance Rate", col="Proportion of complete cases", title=title) 
  
  p + facet_grid(data_type ~ ., switch = "x") +
    theme(panel.spacing = unit(1, "lines")) + # add spacing between plots
    theme(strip.background = element_blank()) +  # remove the background 
    theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5, vjust=1.5)) + 
    theme(axis.title.y = element_text(size=15, vjust=1.75)) + 
    theme(axis.title.x = element_text(size=15)) + 
    theme(legend.title = element_text(size=15)) + 
    theme(legend.text = element_text(size = 14)) +
    theme(legend.position="bottom") +  
    scale_color_manual(labels = c("6%", "30%", "48%"),
                       values = c("aquamarine4", "coral3", "royalblue1")) 
  
}

plot_df_mar_cont = process_df(label="Sequential MNAR does hold", name="SeqMNAR_continuous")
plot_df_per_cont = process_df(label="Sequential MNAR does not hold", name="Permutation_continuous")
plot_d_cont = rbind(plot_df_mar_cont, plot_df_per_cont)

plot_df_mar_bin = process_df(label="Sequential MNAR does hold", name="SeqMNAR_binary")
plot_df_per_bin = process_df(label="Sequential MNAR does not hold", name="Permutation_binary")
plot_d_bin = rbind(plot_df_mar_bin, plot_df_per_bin)

plot_func(plot_d_cont, title="Gaussian Data")
plot_func(plot_d_bin, title="Binary Data")



