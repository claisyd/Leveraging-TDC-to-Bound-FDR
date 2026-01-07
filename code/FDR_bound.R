# FDR_bound_new.R

library(dplyr)

# Helper: s_tau; affine map from (tau, 1/2] to [0, tau] inverted.
# Gives second fold for BD
# Note: p here is first-fold p

s_tau = function(p, tau) {
  
  sp = tau * (1/2 - p) / (1/2 - tau)
  return(sp)
  
}

# p to tuples
# Input: p-values, cutoff tau
# Output: label (T/ED/BD), W1, W2, L1, L2 for EACH p
# |T|, |ED|, |BD| need not be equal but we need n labels and winning scores for index later on

pvalues_to_tuples = function(p, tau) {
  
  n = length(p)
  
  # Labels before the first fold
  # T: [0, tau], BD: (tau, 1 - tau), ED: [1 - tau, 1]
  label = case_when(
    p <= tau            ~ "T",
    p >= 1 - tau        ~ "ED",
    p > tau & p < 1-tau ~ "BD",
    TRUE                ~ NA_character_ 
  )
  
  # Save labels and p-value in one data frame for easier ref
  data = data.frame(
    index = seq_len(n),
    p_value = p,
    label = factor(label, levels = c("T","BD","ED"))
  )
  
  # First fold; find W1 = mirrored p-value
  # so if T or ED -> W1 in [0, tau]
  # if BD -> W1 in (tau, 1/2]
  W1 = pmin(p, 1 - p)
  data = data %>% mutate("W1" = W1)
  
  # Scoring label 1: find L1
  # L1 = 1 if T, -1 if ED, 0 if BD
  L1 = case_when(
    data$label == "T"   ~ 1,
    data$label == "ED"  ~ -1,
    data$label == "BD"  ~ 0,
    TRUE                ~ NA_real_ 
  )
  data = data %>% mutate("L1" = L1)
  
  # Second fold: find W2, mapping from [0, 1/2] to [0, tau].
  # so if T or ED, W2 = W1; else it is using s_tau
  W2 = case_when(
    data$label == "T"   ~ data$W1,
    data$label == "ED"  ~ data$W1,
    data$label == "BD"  ~ s_tau(data$W1, tau),
    TRUE                ~ NA_real_
  )
  data = data %>% mutate("W2" = W2)
  
  # Scoring label 2: Find L2
  # +1 if T or ED (ie, T2), -1 if BD (ie, D2)
  L2 = case_when(
    data$label == "T"   ~ 1,
    data$label == "ED"  ~ 1,
    data$label == "BD"  ~ -1,
    TRUE                ~ NA_real_ 
  )
  data = data %>% mutate("L2" = L2)
  
  # change winning scores to 1 - p
  data = data %>% mutate("W1" = 1 - W1, "W2" = 1 - W2) 
  
  return(data)
}

# FDR_bound
# input: data of (W1, W2, L1, L2) and original labels; alpha, ratio r of E(BD) / E(ED)
# and c = to add onto the omega function
# Output (for now): estFDR, # of discoveries
# Basically should be Eric's code from step 3 onward

FDR_bound = function(data, alpha, r, c = 1) { 
  
  # Sort W2 descending
  sorted_indices = order(data$W2, decreasing = T) 
  sorted_data = data[sorted_indices, ]
  
  L2_sorted = sorted_data$L2
  L1_sorted = sorted_data$L1
  
  # Largest subarray such that omega(m) <= alpha
  T2 = cumsum(L2_sorted == 1)
  D2 = cumsum(L2_sorted == -1)
  
  omega = ((D2 + c) / r) / pmax(1, T2 - (D2/r))
  ok = which(omega <= alpha) 
  
  if (length(ok) == 0L) return(list("FDR_hat" = 0,
                                    "n_targets" = 0,
                                    "n_decoys" = 0,
                                    "selected_indices" = numeric(0),
                                    "cutoff" = NA))
  
  x = max(ok) # cutoff
  
  # Take first x of the sorted data
  T2_x = subset(sorted_data[1:x,], L2 == 1) # pseudotarget in top x sorted W2
  
  # Find n_tar and n_dec from here
  n_tar = sum(T2_x$L1 == 1)
  n_dec = sum(T2_x$L1 == -1)
  
  estFDR = n_dec / max(n_tar, 1)
  
  # Selected indices (for FDP and power)
  targets_selected = T2_x[which(T2_x$L1 == 1), ]
  decoys_selected = T2_x[which(T2_x$L1 == -1), ]
  sel = targets_selected$index
  
  return(list("FDR_hat" = estFDR,
              "n_targets" = n_tar,
              "n_decoys" = n_dec,
              "selected_indices" = sel,
              "cutoff" = x))
  
}

# tdc_from_p
# Wrapper function

tdc_from_p = function(p, alpha, tau, c = 1) {
  
  # p to TDC
  TDC = pvalues_to_tuples(p, tau)
  
  # Ratio of BD to ED 
  r = (1 - 2 * tau) / tau
  
  # FDR_bound
  FDR_bound(TDC, alpha, r, c)
  
}
