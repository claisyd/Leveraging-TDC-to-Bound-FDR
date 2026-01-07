# Plotting helper functions, to trim out the Dependent.qmd code

# Cluster summaries
cluster_summaries = function(df, firstK = 12) {
  
  # p-value density (overlay nulls vs alts)
  p1 = ggplot(df, aes(p, colour = is_null, fill = is_null)) +
    geom_density(alpha = 0.25) +
    scale_color_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#ff7f0e"),
                       labels = c("FALSE" = "Alts", "TRUE" = "Nulls"), name = NULL) +
    scale_fill_manual(values  = c("FALSE" = "#1f77b4", "TRUE" = "#ff7f0e"),
                      labels = c("FALSE" = "Alts", "TRUE" = "Nulls"), name = NULL) +
    labs(title = "p-value density", x = "p", y = "density")
  
  # Boxplots by cluster (first K clusters) 
  df_firstK = df %>% filter(cluster <= firstK)
  
  p2 = ggplot(df_firstK, aes(factor(cluster), p, fill = factor(cluster))) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = "none") +
    labs(title = sprintf("p by cluster (first %d clusters)", firstK),
         x = "cluster", y = "p")
  
  # First 400 z-scores by index (colored by cluster)
  df_head = df %>% slice_head(n = 400)
  p3 = ggplot(df_head, aes(i, z, colour = factor(cluster))) +
    geom_point(size = 1.2, alpha = 0.9) +
    guides(colour = "none") +
    labs(title = "First 400 z-scores by index (colored by cluster)",
         x = "index", y = "z")
  
  # Per cluster summaries
  cluster_stats = df %>%
    group_by(cluster) %>%
    summarise(mean_p = mean(p), alt_rate = 1 - mean(is_null), .groups = "drop")
  
  p4 = ggplot(cluster_stats, aes(alt_rate, mean_p, colour = alt_rate)) +
    geom_point(alpha = 0.7) +
    guides(colour = guide_colorbar(title = "alt rate")) +
    labs(title = "Cluster mean p vs alternative rate",
         x = "alt rate in cluster", y = "mean p")
  
  return(p1 + p2 + p3 + p4)
  
}

# Null histogram plotter
null_histogram = function(n, M, structure, ..., .generator = NULL, title_txt, param_txt) {
  
  hist_breaks = seq(0, 1, length.out = 41) # 40 bins for histogram
  hist_counts = integer(length(hist_breaks) - 1)
  
  for (m in 1:M) {
    sm = generate_p_dependent(n, pi0 = 1, mu_alt = 0, structure, ..., .generator = .generator)
    p  = sm$p  # all nulls here
    
    # histogram aggregation
    hist_counts = hist_counts +
      tabulate(findInterval(p, hist_breaks, rightmost.closed = TRUE),
               nbins = length(hist_breaks) - 1)
    
  }
  
  N_total = M * n
  hist_mids = head(hist_breaks, -1) + diff(hist_breaks)/2
  df_hist = data.frame(p_mid = hist_mids, count = hist_counts)
  
  p_hist = ggplot(df_hist, aes(x = p_mid, y = count)) +
    geom_col(width = diff(hist_breaks)[1], fill = "grey80", color = "grey40") +
    labs(title = title_txt, subtitle = param_txt, x = "p", y = "Frequency") +
    theme_linedraw()
  
  return(p_hist)
  
}

# FDP plot
fdp_plot = function(df, main_txt, params_txt) {
  
  cols = c("FDR_bound" = "#009E73",  # green
           "Storey"    = "#E69F00",  # orange
           "BH"        = "#56B4E9",  # sky blue
           "BY"        = "#CC79A7")  # purple
  
  p_fdp = ggplot(df, aes(x = alpha, y = FDP_mean, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = pmax(0, FDP_qL),
                    ymax = pmin(1, FDP_qU)),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(
      title = main_txt, subtitle = params_txt,
      x = expression(alpha),
      y = "Empirical FDR",
      color = "Method", fill = "Method"
    ) +
    theme_linedraw(base_size = 13) + geom_abline(slope = 1, intercept = 0,
                                                linetype = "dashed", linewidth = 0.6, color = "grey40")
  
  p_fdp
  
}

# Power plot
power_plot = function(df, main_txt, params_txt) {
  
  cols = c("FDR_bound" = "#009E73",  # green
           "Storey"    = "#E69F00",  # orange
           "BH"        = "#56B4E9",  # sky blue
           "BY"        = "#CC79A7")  # purple
  
  p_power = ggplot(df, aes(x = alpha, y = Power_mean, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = pmax(0, Power_qL),
                    ymax = pmin(1, Power_qU)),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(
      title = main_txt, subtitle = params_txt,
      x = expression(alpha),
      y = "Power",
      color = "Method", fill = "Method"
    ) +
    theme_linedraw(base_size = 13)
  
  p_power
  
}

# filter function and adds MCSE
res_filter = function(res) {
  
  df = res %>% filter(Method == "FDR_bound") %>%
    mutate(
      gap = FDR_bound_mean - FDP_mean,
      gap_raw = FDR_bound_raw_mean - FDP_mean,
      
      FDP_mcse_mean  = FDP_sd  / sqrt(M),
      FDRb_mcse_mean = FDR_bound_sd / sqrt(M),
      FDRb_raw_mcse_mean = FDR_bound_raw_sd / sqrt(M),
      gap_mcse = sqrt(FDRb_mcse_mean^2 + FDP_mcse_mean^2),
      gap_raw_mcse = sqrt(FDRb_raw_mcse_mean^2 + FDP_mcse_mean^2)) %>% filter(alpha <= 0.25)
  
}

# difference plot
bound_fdp_plot = function(df, title_txt, params_txt) {
  ggplot(df, aes(x = alpha)) +
    
    # quantile ribbons
    geom_ribbon(aes(ymin = FDP_qL,  ymax = FDP_qU),
                fill = "#009E73", alpha = 0.2, colour = NA, show.legend = FALSE) +
    geom_ribbon(aes(ymin = FDR_bound_qL, ymax = FDR_bound_qU),
               fill = "salmon", alpha = 0.2, colour = NA, show.legend = FALSE) +
    
    # mean lines
    geom_line(aes(y = FDP_mean, linetype = "FDP"),
              linewidth = 1, colour = "#009E73") +
    geom_line(aes(y = FDR_bound_mean, linetype = "FDR Bound"),
              linewidth = 1, colour = "salmon") +
    scale_linetype_manual(values = c("FDR Bound" = "dashed", "FDP" = "solid"),
                          breaks = c("FDP", "FDR Bound"), name = NULL) +
    guides(fill = "none", colour = "none", alpha = "none") +
    
    labs(title = title_txt, subtitle = params_txt,
         x = expression(alpha), y = "Rate") +
    theme_linedraw(base_size = 13) +
    theme(legend.position = "top")
  
}

# Gap plot
gap_plot = function(df, title_txt) {
  
  p1 = ggplot(df, aes(x = alpha, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "goldenrod2") +
    geom_ribbon(aes(ymin = Gap_qL, ymax = Gap_qU),
                fill = "goldenrod2", alpha = 0.2, colour = NA) + 
    labs(
      title = title_txt, subtitle = "Truncated by 1",
      x = expression(alpha), y = "Gap"
    ) +
    theme_linedraw(base_size = 13)
  
  p2 = ggplot(df, aes(x = alpha, y = gap_raw)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "mediumslateblue") +
    geom_ribbon(aes(ymin = Gap_raw_qL, ymax = Gap_raw_qU),
                fill = "mediumslateblue", alpha = 0.2, colour = NA) + 
    labs(
      title = NULL, subtitle = "Raw FDR Bound",
      x = expression(alpha), y = "Gap"
    ) +
    theme_linedraw(base_size = 13)
  
  p1 + p2
  
}

# 1x2 plot for FDP and power

fdp_power_plot = function(df, main_txt) {
  
  p_fdp = ggplot(df, aes(x = alpha, y = FDP_mean, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = pmax(0, FDP_qL),
                    ymax = pmin(1, FDP_qU)),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(
      title = main_txt, subtitle = "FDR vs α",
      x = expression(alpha),
      y = "Empirical FDR",
      color = "Method", fill = "Method"
    ) +
    theme_linedraw(base_size = 13) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
    theme(legend.position = "none")
  
  p_power =   p_power = ggplot(df, aes(x = alpha, y = Power_mean, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = pmax(0, Power_qL),
                    ymax = pmin(1, Power_qU)),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(
      title = NULL, subtitle = "Power vs α",
      x = expression(alpha),
      y = "Power",
      color = "Method", fill = "Method"
    ) +
    theme_linedraw(base_size = 13)
  
  p_fdp + p_power
  
}

# Gap plot 
gap_plot_single = function(df, title_txt) {
  
  p1 = ggplot(df, aes(x = alpha, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "goldenrod2") +
    geom_ribbon(aes(ymin = Gap_qL, ymax = Gap_qU),
                fill = "goldenrod2", alpha = 0.2, colour = NA) + 
    labs(
      title = title_txt,
      x = expression(alpha), y = "Difference"
    ) +
    theme_linedraw(base_size = 13)
  
  p1
  
}

# Gap plot without bands
gap_plot_no_bands = function(df, title_txt) {
  
  p1 = ggplot(df, aes(x = alpha, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "goldenrod2") +
    labs(
      title = title_txt,
      x = expression(alpha), y = "Mean FDR Bound - Mean FDP"
    ) +
    theme_linedraw(base_size = 13)
  
  p1
  
}

# Gap plot
gap_plot2_no_bands = function(df, title_txt) {
  
  p1 = ggplot(df, aes(x = alpha, y = gap)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "goldenrod2") +
    labs(
      title = title_txt, subtitle = "Clipped at 1",
      x = expression(alpha), y = "Gap"
    ) +
    theme_linedraw(base_size = 13)
  
  p2 = ggplot(df, aes(x = alpha, y = gap_raw)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1, color = "mediumslateblue") +
    labs(
      title = NULL, subtitle = "Raw FDR Bound",
      x = expression(alpha), y = "Gap"
    ) +
    theme_linedraw(base_size = 13)
  
  p1 + p2
  
}