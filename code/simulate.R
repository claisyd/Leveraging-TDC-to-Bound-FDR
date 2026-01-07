# Simulation helper functions

library(dplyr)
library(knitr)

# generate_p
generate_p = function(n, p, u_min = 0, u_max = 1,
                      beta_a = 0.5, beta_b = 5, clamp_to_open01 = TRUE) {
  
  stopifnot(length(n) == 1, n >= 1,
            length(p) == 1, is.numeric(p), p >= 0, p <= 1,
            is.numeric(u_min), is.numeric(u_max), u_min < u_max,
            u_min >= 0, u_max <= 1,
            is.numeric(beta_a), is.numeric(beta_b), beta_a > 0, beta_b > 0)
  
  z  = rbinom(n, 1, p)            # 1 => Uniform (TRUE NULL), 0 => Beta (FALSE NULL)
  is_null = (z == 1L)
  x  = numeric(n)
  kU = sum(is_null)
  
  if (kU > 0L) x[is_null]  = runif(kU, min = u_min, max = u_max)
  if (kU < n)  x[!is_null] = rbeta(n - kU, shape1 = beta_a, shape2 = beta_b)
  
  if (clamp_to_open01) {
    eps = .Machine$double.eps
    x = pmin(pmax(x, eps), 1 - eps)  # safer for qvalue later
  }
  
  list(
    p       = x,
    is_null = is_null,
    is_alt  = !is_null,
    labels  = factor(ifelse(is_null, "null", "alt"), levels = c("null","alt"))
  )
}

# iterate
iterate = function(n, alpha, tau, params,
                    structure = c("independent", "weak_local", "strong_local", "null_clumps", "anti_pairs"),
                    dep_mu_alt = 2, anti_rho = 0.85, c = 1,
                    block_size = 10,
                    .generator = NULL      # optional custom dependent generator
) {
  structure = match.arg(structure)
  
  w     = params[1]
  u_min = params[2]
  u_max = params[3]
  a     = params[4]
  b     = params[5]
  
  if (structure == "independent") {
    # your existing independent generator (keep as-is)
    sim = generate_p(
      n = n, p = w,
      u_min = u_min, u_max = u_max,
      beta_a = a, beta_b = b,
      clamp_to_open01 = TRUE
    )
  } else {
    # refactored dependent generators via wrapper / custom function
    sim = generate_p_dependent(
      n = n, pi0 = w, mu_alt = dep_mu_alt,
      structure = structure,
      anti_rho = anti_rho, block_size = block_size,
      .generator = .generator
    )
  }
  
  p       = sim$p
  is_null = sim$is_null
  
  # Methods
  fb = metrics_FDR_bound(p, alpha, tau, is_null, c)   # must provide $FDR_bound_raw (and optionally $FDR_bound)
  bh = metrics_BH(p, is_null, alpha)
  by = metrics_BY(p, is_null, alpha)
  st = metrics_Storey(p, is_null, alpha)
  
  FDRb_raw = fb$FDR_bound_raw
  FDRb     = if (!is.null(fb$FDR_bound)) fb$FDR_bound else pmin(FDRb_raw, 1)
  
  list(
    # FDR-bound procedure outputs
    FB_FDP   = fb$FDP,
    FB_Power = fb$Power,
    FB_R     = fb$n_targets,
    FB_FDRb_raw = FDRb_raw,
    FB_FDRb     = FDRb,
    
    # Comparators
    BH_FDP   = bh$FDP,   BH_Power = bh$Power,   BH_R = bh$R,
    BY_FDP   = by$FDP,   BY_Power = by$Power,   BY_R = by$R,
    ST_FDP   = st$FDP,   ST_Power = st$Power,   ST_R = st$R
  )
}

# simulate
simulate = function(n, alpha, tau, params, B,
                     structure = c("independent", "weak_local", "strong_local", "null_clumps", "anti_pairs"),
                     dep_mu_alt = 2, anti_rho = 0.85, diagn = FALSE, c = 1,
                     block_size = 10,
                     .generator = NULL,      # optional custom dependent generator
                    qL = 0.25, qU = 0.75
) {
  structure = match.arg(structure)
  
  # storage
  res_FB_FDP = res_FB_Power = res_FB_R = numeric(B)
  res_FB_FDRb_raw = res_FB_FDRb = numeric(B)
  res_BH_FDP = res_BH_Power = res_BH_R = numeric(B)
  res_BY_FDP = res_BY_Power = res_BY_R = numeric(B)
  res_ST_FDP = res_ST_Power = res_ST_R = numeric(B)
  
  for (i in seq_len(B)) {
    out = iterate(
      n, alpha = alpha, tau = tau, params = params,
      structure = structure, dep_mu_alt = dep_mu_alt, c = c, anti_rho = anti_rho,
      block_size = block_size, .generator = .generator
    )
    # FDR-bound
    res_FB_FDP[i]      = out$FB_FDP
    res_FB_Power[i]    = out$FB_Power
    res_FB_R[i]        = out$FB_R
    res_FB_FDRb_raw[i] = out$FB_FDRb_raw
    res_FB_FDRb[i]     = out$FB_FDRb
    # BH
    res_BH_FDP[i] = out$BH_FDP; res_BH_Power[i] = out$BH_Power; res_BH_R[i] = out$BH_R
    # BY
    res_BY_FDP[i] = out$BY_FDP; res_BY_Power[i] = out$BY_Power; res_BY_R[i] = out$BY_R
    # Storey
    res_ST_FDP[i] = out$ST_FDP; res_ST_Power[i] = out$ST_Power; res_ST_R[i] = out$ST_R
  }
  
  # empirical quantiles across replicates (for FDP of all methods, and power)
  # Our method's FDP and power
  FB_qL   = as.numeric(quantile(res_FB_FDP, qL, na.rm = TRUE))
  FB_qU   = as.numeric(quantile(res_FB_FDP, qU, na.rm = TRUE))
  FB_power_qL = as.numeric(quantile(res_FB_Power, qL, na.rm = TRUE))
  FB_power_qU = as.numeric(quantile(res_FB_Power, qU, na.rm = TRUE))
  
  # Bound 
  Bound_qL  = as.numeric(quantile(res_FB_FDRb, qL, na.rm = TRUE))
  Bound_qU  = as.numeric(quantile(res_FB_FDRb, qU, na.rm = TRUE))
  
  # BH and power
  BH_qL  = as.numeric(quantile(res_BH_FDP, qL, na.rm = TRUE))
  BH_qU  = as.numeric(quantile(res_BH_FDP, qU, na.rm = TRUE))
  BH_power_qL = as.numeric(quantile(res_BH_Power, qL, na.rm = TRUE))
  BH_power_qU = as.numeric(quantile(res_BH_Power, qU, na.rm = TRUE))
  
  # BY and power
  BY_qL  = as.numeric(quantile(res_BY_FDP, qL, na.rm = TRUE))
  BY_qU  = as.numeric(quantile(res_BY_FDP, qU, na.rm = TRUE))
  BY_power_qL = as.numeric(quantile(res_BY_Power, qL, na.rm = TRUE))
  BY_power_qU = as.numeric(quantile(res_BY_Power, qU, na.rm = TRUE))
  
  # Storey and power
  ST_qL  = as.numeric(quantile(res_ST_FDP, qL, na.rm = TRUE))
  ST_qU  = as.numeric(quantile(res_ST_FDP, qU, na.rm = TRUE))
  ST_power_qL = as.numeric(quantile(res_ST_Power, qL, na.rm = TRUE))
  ST_power_qU = as.numeric(quantile(res_ST_Power, qU, na.rm = TRUE))
  
  # Gap
  gap = res_FB_FDRb - res_FB_FDP
  gap_raw = res_FB_FDRb_raw - res_FB_FDP
  gap_qL = as.numeric(quantile(gap, qL, na.rm = TRUE))
  gap_qU = as.numeric(quantile(gap, qU, na.rm = TRUE))
  gap_raw_qL = as.numeric(quantile(gap_raw, qL, na.rm = TRUE))
  gap_raw_qU = as.numeric(quantile(gap_raw, qU, na.rm = TRUE))
  
  # Table
  data.frame(
    Method   = c("FDR_bound", "BH", "BY", "Storey"),
    FDP_mean = c(mean(res_FB_FDP), mean(res_BH_FDP), mean(res_BY_FDP), mean(res_ST_FDP, na.rm = TRUE)),
    FDP_sd   = c(sd(res_FB_FDP),   sd(res_BH_FDP),   sd(res_BY_FDP),   sd(res_ST_FDP, na.rm = TRUE)),
    Power_mean = c(mean(res_FB_Power), mean(res_BH_Power), mean(res_BY_Power), mean(res_ST_Power, na.rm = TRUE)),
    Power_sd   = c(sd(res_FB_Power),   sd(res_BH_Power),   sd(res_BY_Power),   sd(res_ST_Power, na.rm = TRUE)),
    R_mean     = c(mean(res_FB_R),     mean(res_BH_R),     mean(res_BY_R),     mean(res_ST_R, na.rm = TRUE)),
    R_sd       = c(sd(res_FB_R),       sd(res_BH_R),       sd(res_BY_R),       sd(res_ST_R, na.rm = TRUE)),
    
    # Empirical quantile ribbons for all methods
    FDP_qL = c(FB_qL, BH_qL, BY_qL, ST_qL),
    FDP_qU = c(FB_qU, BH_qU, BY_qU, ST_qU),
    
    Power_qL = c(FB_power_qL, BH_power_qL, BY_power_qL, ST_power_qL),
    Power_qU = c(FB_power_qU, BH_power_qU, BY_power_qU, ST_power_qU),
    
    # Bound summaries on FB row
    FDR_bound_mean     = c(mean(res_FB_FDRb), NA, NA, NA),
    FDR_bound_sd       = c(sd(res_FB_FDRb),   NA, NA, NA),
    FDR_bound_raw_mean = c(mean(res_FB_FDRb_raw), NA, NA, NA),
    FDR_bound_raw_sd   = c(sd(res_FB_FDRb_raw),   NA, NA, NA),
    Gap_bound_minus_FDP = c(mean(res_FB_FDRb) - mean(res_FB_FDP), NA, NA, NA),
    
    FDR_bound_qL = c(Bound_qL, NA, NA, NA),
    FDR_bound_qU = c(Bound_qU, NA, NA, NA),
    
    Gap_qL = c(gap_qL, NA, NA, NA),
    Gap_qU = c(gap_qU, NA, NA, NA),
    Gap_raw_qL = c(gap_raw_qL, NA, NA, NA),
    Gap_raw_qU = c(gap_raw_qU, NA, NA, NA),
    
    Structure = structure,
    alpha     = alpha,
    n         = n,
    pi0       = params[1],
    tau       = tau,
    c         = c,
    iters     = B,
    Storey_NA = c("-", "-", "-", sum(is.na(res_ST_FDP)))
  )
}

# summary table
summary_table = function(summary_df, caption = "Simulation summary", digits = 3) {
  stopifnot(is.data.frame(summary_df))
  
  # Desired columns in this order with new display names
  col_map = c(
    Method     = "Method",
    FDP_mean   = "FDP",
    FDP_sd     = "FDP SD",
    Power_mean = "Power",
    Power_sd   = "Power SD",
    R_mean     = "Discoveries",
    R_sd       = "Discoveries SD",
    Storey_NA  = "NA Count"
  )
  
  # Keep only columns that exist; build display data.frame
  keep = intersect(names(col_map), names(summary_df))
  if (length(keep) == 0L) stop("No expected columns found in `summary_df`.")
  
  out = summary_df[keep]
  
  # Round numeric columns (leave character/factor as-is)
  is_num = vapply(out, is.numeric, logical(1))
  out[is_num] = lapply(out[is_num], function(x) round(x, digits))
  
  # Rename columns for display
  names(out) = unname(col_map[keep])
  
  knitr::kable(out, caption = caption, booktabs = TRUE, align = "lccccccc")
}