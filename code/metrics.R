# Metrics for BH / BY / Storey

library(qvalue)

# BH procedure and metrics
metrics_BH = function(p, is_null, alpha) {
  
  stopifnot(length(p) == length(is_null), is.logical(is_null))
  adj = p.adjust(p, method = "BH")
  rejects = (adj <= alpha)
  
  R  = sum(rejects)
  FD = if (R == 0L) 0L else sum(rejects & is_null)
  TP = R - FD
  m1 = sum(!is_null)
  
  FDP   = if (R == 0L) 0 else FD / R
  Power = if (m1 == 0L) NA_real_ else TP / m1
  est_FDR_set = if (R == 0L) 0 else max(adj[rejects])
  
  list(R = R, FD = FD, TP = TP, FDP = FDP, Power = Power,
       rejects = rejects, est_FDR_set = est_FDR_set, adj_p = adj)
  
}

# BY procedure and metrics
metrics_BY = function(p, is_null, alpha) {
  
  stopifnot(length(p) == length(is_null), is.logical(is_null))
  adj = p.adjust(p, method = "BY")
  rejects = (adj <= alpha)
  
  R  = sum(rejects)
  FD = if (R == 0L) 0L else sum(rejects & is_null)
  TP = R - FD
  m1 = sum(!is_null)
  
  FDP   = if (R == 0L) 0 else FD / R
  Power = if (m1 == 0L) NA_real_ else TP / m1
  est_FDR_set = if (R == 0L) 0 else max(adj[rejects])
  
  list(R = R, FD = FD, TP = TP, FDP = FDP, Power = Power,
       rejects = rejects, est_FDR_set = est_FDR_set, adj_p = adj)
  
}

# Storey procedure and metrics (fail-safe NA on any error)
metrics_Storey = function(p, is_null, alpha) {
  tryCatch({
    # basic checks (kept inside tryCatch so we still return NAs if they fail)
    stopifnot(length(p) == length(is_null), is.logical(is_null), length(alpha) == 1L)
    
    # clamp p into (0,1) to avoid boundary weirdness
    eps = .Machine$double.eps
    p2  = pmin(pmax(p, eps), 1 - eps)
    
    # plain-vanilla Storey q-values — no parameter tweaks, no fallback
    qv = qvalue(p2)$qvalues
    
    # if qvalue returned any NAs, treat as failure
    if (anyNA(qv)) stop("qvalue returned NA q-values")
    
    rejects = (qv <= alpha)
    
    R  = sum(rejects)
    FD = if (R == 0L) 0L else sum(rejects & is_null)
    TP = R - FD
    m1 = sum(!is_null)
    
    FDP         = if (R == 0L) 0 else FD / R
    Power       = if (m1 == 0L) NA_real_ else TP / m1
    est_FDR_set = if (R == 0L) 0 else max(qv[rejects])
    
    list(
      R = R, FD = FD, TP = TP,
      FDP = FDP, Power = Power,
      rejects = rejects,
      est_FDR_set = est_FDR_set,
      qvalues = qv
    )
  },
  error = function(e) {
    # On any failure, return NAs of the right types/shapes
    list(
      R = NA_integer_, FD = NA_integer_, TP = NA_integer_,
      FDP = NA_real_, Power = NA_real_,
      rejects = rep(NA, length(p)),
      est_FDR_set = NA_real_,
      qvalues = rep(NA_real_, length(p))
    )
  })
}

# FDR Bound procedure and metrics
metrics_FDR_bound = function(p, alpha, tau, is_null, c = 1) {
  
  TDC = pvalues_to_tuples(p, tau) %>% mutate("null" = is_null)
  r = (1 - 2 * tau) / tau
  out = FDR_bound(TDC, alpha, r, c)
  sel = out$selected_indices
  
  nT = out$n_targets # Discoveries/Targets
  FD  = if (nT == 0L) 0L else sum(is_null[sel]) # False discoveries
  TP  = nT - FD # True positives
  m1  = sum(!is_null) # number of alts
  
  FDP   = if (nT == 0L) 0 else FD / nT
  Power = if (m1 == 0L) NA_real_ else TP / m1
  
  FDR_raw = out$FDR_hat
  FDR = min(1, FDR_raw) # truncate
  
  return(list(
    "FDP" = FDP,
    "Power" = Power,
    "n_targets" = nT, 
    "FDR_bound" = FDR,
    "FDR_bound_raw" = FDR_raw))
  
}
