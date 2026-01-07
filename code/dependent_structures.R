# Dependent data generation

library(mvtnorm)

# Structure 1-4: Anti-pairs, Null Clumps, Weak local, and strong local
# Helper: Zhang et al. weak-dependence 10x10 covariance (ρ=±0.25)
.cov_zhang10 = function() {
  S = diag(10)
  # first 5: +0.25 among themselves
  S[1:5, 1:5] = 0.25
  diag(S[1:5, 1:5]) = 1
  # cross (first 5) x (last 5): -0.25
  S[1:5, 6:10] = -0.25
  S[6:10, 1:5] = -0.25
  # last 5 among themselves: 0 off-diagonal (already set)
  S
}


# Helper: common return builder
.make_return = function(p, is_null, z, cluster_id, structure, extras = list()) {
  list(
    p = p,
    is_null = is_null,
    z = z,
    cluster_id = cluster_id,
    params_used = c(list(structure = structure), extras)
  )
}

# 1) weak_local: Zhang-style 10x10 covariance blocks (needs mvtnorm::rmvnorm and your .cov_zhang10())
gen_weak_local = function(n, pi0 = 0.8, mu_alt = 2, ...) {

  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  z = numeric(n)
  
  S10 = .cov_zhang10()
  nb = n %/% 10
  rem = n %% 10
  
  if (nb > 0) {
    Zmat = mvtnorm::rmvnorm(nb, sigma = S10)   # nb x 10
    z[seq_len(nb * 10)] = as.vector(t(Zmat))
  }
  if (rem > 0) z[(nb * 10 + 1):n] = rnorm(rem)
  
  z[is_alt] = z[is_alt] + mu_alt
  p = pnorm(z, lower.tail = FALSE)
  cluster_id = rep(seq_len(ceiling(n / 10)), each = 10)[1:n]
  
  .make_return(p, is_null, z, cluster_id, "weak_local",
               extras = list(pi0 = pi0, mu_alt = mu_alt, block_size = 10))
}

# 2) strong_local: blocks of 5 with identical z within each block
gen_strong_local = function(n, pi0 = 0.8, mu_alt = 2, ...) {

  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  z = numeric(n)
  
  nb = n %/% 5
  rem = n %% 5
  if (nb > 0) {
    basez = rnorm(nb)
    Zmat  = matrix(rep(basez, each = 5), nrow = nb, ncol = 5, byrow = TRUE)
    z[seq_len(nb * 5)] = as.vector(t(Zmat))
  }
  if (rem > 0) z[(nb * 5 + 1):n] = rnorm(rem)
  
  z[is_alt] = z[is_alt] + mu_alt
  p = pnorm(z, lower.tail = FALSE)
  cluster_id = rep(seq_len(ceiling(n / 5)), each = 5)[1:n]
  
  .make_return(p, is_null, z, cluster_id, "strong_local",
               extras = list(pi0 = pi0, mu_alt = mu_alt, block_size = 5))
}

# 3) null_clumps: comonotone blocks among nulls (stress Storey’s pi0-hat)
gen_null_clumps = function(n, pi0 = 0.8, mu_alt = 2, block_size = 10, ...) {
  
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  z = rnorm(n)  # start independent
  
  null_idx = which(is_null)
  m0 = length(null_idx)
  if (m0 > 0) {
    B = max(1, ceiling(m0 / block_size))
    blocks = split(null_idx, rep(1:B, length.out = m0))
    for (b in seq_along(blocks)) {
      Ub = runif(1)         # shared null p-value
      Zb = qnorm(1 - Ub)    # so that p = 1 - Phi(Zb) = Ub
      z[blocks[[b]]] = Zb
    }
  }
  
  z[is_alt] = z[is_alt] + mu_alt
  p = pnorm(z, lower.tail = FALSE)
  
  # cluster ids only for nulls' clumps
  cluster_id = rep(NA_integer_, n)
  if (m0 > 0) {
    B = max(1, ceiling(m0 / block_size))
    cluster_id[null_idx] = rep(seq_len(B), length.out = m0)
  }
  
  .make_return(p, is_null, z, cluster_id, "null_clumps",
               extras = list(pi0 = pi0, mu_alt = mu_alt, block_size = block_size))
}

# 4) anti_pairs: null–null pairs with Corr = -anti_rho (others i.i.d.)
gen_anti_pairs = function(n, pi0 = 0.8, mu_alt = 2, anti_rho = 0.85, ...) {
  
  stopifnot(anti_rho >= 0, anti_rho <= 1)
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  eps_all = rnorm(n)
  z = eps_all
  
  K = n %/% 2
  if (anti_rho > 0 && K > 0) {
    s1 = sqrt(anti_rho); s2 = sqrt(1 - anti_rho)
    G  = rnorm(K)
    for (k in seq_len(K)) {
      i = 2 * k - 1; j = i + 1
      if (is_null[i] && is_null[j]) {
        z[i] = s1 * G[k] + s2 * eps_all[i]
        z[j] = -s1 * G[k] + s2 * eps_all[j]
      }
    }
  }
  
  z[is_alt] = z[is_alt] + mu_alt
  p = pnorm(z, lower.tail = FALSE)
  cluster_id = rep(seq_len(ceiling(n / 2)), each = 2)[1:n]
  
  .make_return(p, is_null, z, cluster_id, "anti_pairs",
               extras = list(pi0 = pi0, mu_alt = mu_alt, anti_rho = anti_rho))
}

## ========= Registry + wrapper (backward compatible) =========

GEN_REGISTRY = list(
  weak_local   = gen_weak_local,
  strong_local = gen_strong_local,
  null_clumps  = gen_null_clumps,
  anti_pairs   = gen_anti_pairs
)

# Backward-compatible wrapper: call by structure OR pass a custom .generator function
generate_p_dependent = function(
    n,
    pi0 = 0.8,
    mu_alt = 2,
    structure = c("weak_local", "strong_local", "anti_pairs", "null_clumps"),
    ...,
    .generator = NULL   # optional: a function(n, pi0, mu_alt, ...) returning list(p,is_null, z, cluster_id, params_used)
) {
  if (!is.null(.generator)) {
    return(.generator(n = n, pi0 = pi0, mu_alt = mu_alt, ...))
  }
  structure = match.arg(structure)
  gen = GEN_REGISTRY[[structure]]
  if (is.null(gen)) stop("Unknown structure: ", structure)
  gen(n = n, pi0 = pi0, mu_alt = mu_alt, ...)
}

# Helper: convert generator output -> tidy tibble with cluster + position
to_tibble = function(g, structure) {
  mu = g$params_used$mu_alt
  n  = length(g$p)
  blocksize = if (structure == "weak_local") 10L else if (structure == "strong_local") 5L else 2L
  
  tibble(
    i         = seq_len(n),
    cluster   = g$cluster_id,
    pos_in_blk= ((i - 1L) %% blocksize) + 1L,
    z         = g$z,
    # remove the shift from alts so z_base reflects the *dependence layer*
    z_base    = ifelse(g$is_null, g$z, g$z - mu),
    p         = g$p,
    is_null   = g$is_null,
    structure = structure,
    blocksize = blocksize
  )
}

# Structure 4: Mixed Factor
gen_mixfactor = function(n, pi0 = 0.8, mu_alt, gamma, ...) {
  
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  FF = rnorm(1)
  sgn = sample(c(-1, 1), n, replace = TRUE)
  
  eps = rnorm(n, sd = sqrt(1 - gamma^2))
  z = sgn * gamma * FF + eps 
  
  if (any(is_alt)) z[is_alt] = z[is_alt] + mu_alt
  
  p = pnorm(z, lower.tail = FALSE)
  list(p = p, is_null = is_null, z = z)
}

# Structure 5: Mirror Nulls
gen_mirror_nulls = function(
    n, pi0, mu_alt,
    a = 0.4, b = 4, alt_type = c("beta","normal"), ...
) {
  alt_type = match.arg(alt_type)
  
  p = numeric(n)
  is_null = rep(FALSE, n)
  idx = sample.int(n)        # shuffle roles
  
  m0 = round(pi0 * n)
  
  # All-null case: return early (no alts)
  if (m0 >= n) {
    # pair as many as possible
    m0_pairs = 2 * floor(n / 2)
    null_idx = idx[seq_len(m0_pairs)]
    
    # mirror pairs
    if (m0_pairs > 0) {
      for (k in seq(1, m0_pairs, by = 2)) {
        U = runif(1)
        i = null_idx[k]; j = null_idx[k + 1]
        p[i] = U; p[j] = 1 - U
        is_null[c(i, j)] = TRUE
      }
    }
    # leftover single (if n is odd): give it a Uniform null p
    if (n %% 2 == 1) {
      lone = idx[n]
      p[lone] = runif(1)
      is_null[lone] = TRUE
    }
    
    return(list(
      p = p,
      is_null = is_null,
      z = qnorm(1 - p),
      cluster_id = NULL,
      params_used = list(structure = "mirror_nulls", pi0 = 1, a = a, b = b, alt_type = alt_type)
    ))
  }
  
  # Mixed null/alt case
  # ensure an even number of nulls for pairing
  if (m0 %% 2 == 1) m0 = m0 - 1
  null_idx = if (m0 > 0) idx[seq_len(m0)] else integer(0)
  
  # mirror pairs for nulls
  if (m0 > 0) {
    for (k in seq(1, m0, by = 2)) {
      U = runif(1)
      i = null_idx[k]; j = null_idx[k + 1]
      p[i] = U; p[j] = 1 - U
      is_null[c(i, j)] = TRUE
    }
  }
  
  # alternatives (only if any)
  if (m0 < n) {
    alt_idx = idx[(m0 + 1):n]  # safe because m0 < n here
    if (length(alt_idx) > 0) {
      if (alt_type == "beta") {
        b_alt = rbeta(length(alt_idx), a, b)
        p[alt_idx] = b_alt
      } else {
        z = rnorm(length(alt_idx), mean = mu_alt, sd = 1)
        p[alt_idx] = pnorm(z, lower.tail = FALSE)
      }
    }
  }
  
  list(
    p = p,
    is_null = is_null,
    z = qnorm(1 - p),
    cluster_id = NULL,
    params_used = list(structure = "mirror_nulls", pi0 = pi0, a = a, b = b, alt_type = alt_type)
  )
}

# Structure 6: Clayton Copula
gen_clayton = function(n, pi0 = 0.8, mu_alt = 2.3,
                       theta, alt_type = c("normal","beta"),
                       a = 0.4, b = 4, ...) {
  
  alt_type = match.arg(alt_type)
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  p = vector(mode = "double", length = n)
  
  # Marshall–Olkin sampler for a single Clayton copula over all nulls
  m0 = sum(is_null)
  if (m0 > 0) {
    S  = rgamma(1, shape = 1/theta, rate = 1)  # shape = 1/theta, scale = 1
    Ei = rexp(m0)
    U  = (1 + Ei / S)^(-1/theta)
    p[is_null] = U
  }
  
  # alternatives
  if (any(is_alt)) {
    if (alt_type == "normal") {
      z_alt = rnorm(sum(is_alt), mean = mu_alt, sd = 1)
      p[is_alt] = pnorm(z_alt, lower.tail=FALSE)
    } else {
      b_alt = rbeta(sum(is_alt), a, b)
      p[is_alt] = b_alt
    }
  }
  
  list(
    p = p,
    is_null = is_null,
    z = qnorm(p, lower.tail = FALSE),
    cluster_id = rep(1L, n),
    params_used = list(structure = "clayton", pi0 = pi0, theta = theta, mu_alt = mu_alt)
  )
}

# Structure 7: One Null
gen_one_null = function(n, pi0, a, b) {
  
  # labels
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  # initial vector
  p = numeric(n)
  
  # Set p-values
  null_p = runif(1, 0, 1)
  if (sum(is_null) > 0) {p[is_null] = null_p}
  if (sum(is_alt) > 0) {p[is_alt] = rbeta(sum(is_alt), a, b)}
  
  list(
    p = p,
    is_null = is_null,
    cluster_id = rep(1L, n),
    params_used = list(structure = "one_null", pi0 = pi0, a = a, b = b)
  )
  
}

# Structure 8: One Factor (Dependent Alts)
gen_onefactor = function(n, pi0 = 0.8, mu_alt, gamma, ...) {
  
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  FF = rnorm(1)
  eps = rnorm(n, sd = sqrt(1 - gamma^2))
  z = gamma * FF + eps
  
  if (any(is_alt)) z[is_alt] = z[is_alt] + mu_alt
  
  p = pnorm(z, lower.tail = FALSE)
  list(p = p, is_null = is_null, z = z)
  
}

# Structure 9: One Factor (with Independent Alts)
gen_onefactor_ind_alts = function(n, pi0 = 0.8, mu_alt, gamma, ...) {
  
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  
  FF = rnorm(1)
  eps = rnorm(n, sd = sqrt(1 - gamma^2))
  z = gamma * FF + eps
  
  if (any(is_alt)) z[is_alt] = rnorm(sum(is_alt), mu_alt, 1)
  
  p = pnorm(z, lower.tail = FALSE)
  list(p = p, is_null = is_null, z = z)
  
}

# Structure 10: Block equi-correlated negative Gaussian
gen_neg_equi_block = function(n, pi0 = 0.8, mu_alt, block_size,
                              rho = NULL, ...) {
  
  stopifnot(n >= 1, pi0 >= 0, pi0 <= 1, block_size >= 2)
  if (is.null(rho)) rho = -1 / (block_size - 1)  # strongest feasible negative equicorr in block
  # avoid singular Sigma at the boundary
  if (rho <= -1 / (block_size - 1)) rho = -1 / (block_size - 1) + 1e-6
  stopifnot(rho < 0, rho > -1)  # basic sanity
  
  is_null = rbinom(n, 1, pi0) == 1
  is_alt  = !is_null
  z = numeric(n)
  
  # Generate NULL z-scores: block-wise negative equicorrelation, leftovers IID N(0,1)
  null_idx = which(is_null)
  m0 = length(null_idx)
  if (m0 > 0) {
    # precompute full-block covariance
    J = matrix(1, nrow = block_size, ncol = block_size)
    Sigma_b = (1 - rho) * diag(block_size) + rho * J  # diag=1, offdiag=rho
    
    k_full = m0 %/% block_size
    r_tail = m0 %% block_size
    pos = 1
    
    # full blocks
    if (k_full > 0) {
      for (k in 1:k_full) {
        z_block = as.numeric(rmvnorm(1, mean = rep(0, block_size), sigma = Sigma_b))
        idx_block = null_idx[pos:(pos + block_size - 1)]
        z[idx_block] = z_block
        pos = pos + block_size
      }
    }
    # leftover nulls: IID N(0,1)
    if (r_tail > 0) {
      idx_tail = null_idx[pos:(pos + r_tail - 1)]
      z[idx_tail] = rnorm(r_tail, mean = 0, sd = 1)
    }
  }
  
  # Generate ALT z-scores: IID N(mu_alt, 1), independent of nulls 
  alt_idx = which(is_alt)
  m1 = length(alt_idx)
  if (m1 > 0) {
    z[alt_idx] = rnorm(m1, mean = mu_alt, sd = 1)
  }
  
  p = pnorm(z, lower.tail = FALSE)
  list(p = p, is_null = is_null, z = z)
}

# Structure 11: Guo-Rao Worst Case (Independent Alts)
gen_guo_rao_worstcase = function(n, pi0, alpha, ...) {
  
  stopifnot(n >= 1, pi0 >= 0, pi0 <= 1, alpha > 0, alpha < 1)
  
  m  = n
  m0 = floor(pi0 * m)
  m1 = m - m0
  
  # Randomly choose which indices are true nulls
  # Permutes the indices, and sets the first 800 permuted ones as null
  perm = sample.int(m)
  is_null = rep(FALSE, m)
  if (m0 > 0) is_null[perm[1:m0]] = TRUE
  I0 = which(is_null)
  I1 = which(!is_null)
  
  # Build the U_1,...,U_m, U_{m+1} grid uniforms (width alpha/m)
  U = numeric(m + 1)
  if (m > 0) {
    U[1:m]   = (0:(m-1)) * (alpha/m) + runif(m) * (alpha/m)
  }
  U[m + 1] = runif(1, alpha, 1)
  
  # Probabilities for N (eq. (38) in Guo & Rao)
  H_m0  = if (m0 > 0) sum(1 / (1:m0)) else 0 # Harmonic sum
  alpha_max = 1 / ( (m1/m) + (m0/m)*H_m0 ) # only if p_rest is negative
  probs = numeric(m + 1)
  if (m0 > 0) probs[1:m0] = (m0 / m) * alpha * (1 / (1:m0))      # n = 1..m0, depends on n
  if (m0 + 1 <= m) probs[(m0 + 1):m] = alpha / m                 # n = m0+1..m, fixed
  p_rest = 1 - alpha * (m1 / m + (m0 / m) * H_m0)                # n = m+1, fixed
  if (p_rest < 0) {p_rest = 1 - alpha_max * (m1 / m + (m0 / m) * H_m0)}
  probs[m + 1] = p_rest
  # Note all probs sum = 1 as required
  
  # Draw N and assign p-values per the three phases
  N = sample.int(m + 1, size = 1, prob = probs) # draw a single N for three phases below
  p = rep(NA_real_, m)
  
  if (N <= m0 && m0 > 0) {
    # Phase 1: pick N true nulls to equal U_N; everyone else gets U_{m+1}
    chosen = sample(I0, N, replace = FALSE)
    p[chosen]                 = U[N]
    p[setdiff(1:m, chosen)]   = U[m + 1]
  } else if (N >= m0 + 1 && N <= m) {
    # Phase 2: all m0 nulls + (m - N) alts equal U_N; remaining alts get U_{m+1}
    k_alt = m - N
    chosen_alts = if (k_alt > 0) sample(I1, k_alt, replace = FALSE) else integer(0)
    chosen = c(I0, chosen_alts)
    p[chosen]                 = U[N]
    p[setdiff(1:m, chosen)]   = U[m + 1]
  } else {
    # Phase 3: all nulls get U_{m+1} (>= alpha); all alts get U_1 (<= alpha/m)
    if (length(I0)) p[I0] = U[m + 1]
    if (length(I1)) p[I1] = U[1]
  }
  
  is_alt = !is_null
  p[is_alt] = rbeta(sum(is_alt), 0.4, 4)
  
  list(p = p, is_null = is_null)
}

# Structure 12: Guo-Rao Worst-Case v2 (Alts are zero)
harmonic = function(n) {
  if (n < 1) return(0)
  sum(1 / seq_len(n))
}

.to_probs = function(m0, m1, alpha) {
  # Same mixing law as Guo–Rao (ensures true-null marginals are Uniform)
  m = m0 + m1
  p = numeric(m + 1L)
  if (m0 > 0L) {
    p[1:m0] = (m0 / m) * alpha * (1 / seq_len(m0))
  }
  if (m0 < m) {
    p[(m0 + 1L):m] = alpha / m
  }
  p[m + 1L] = 1 - sum(p[1:m])
  if (p[m + 1L] < -1e-12) {
    stop(glue("alpha of {alpha} too large: leftover probability {p[m+1L]} < 0"))
  }
  p[m + 1L] = max(p[m + 1L], 0)
  p
}

# One draw where only the true-null block is dependent;
# the false-null block is *independent of it* (here: deterministic zeros).
# Returns: list(p, N, is_true, U)
sample_true_only = function(m0, m1, alpha) {
  stopifnot(m0 >= 0L, m1 >= 0L, alpha > 0, alpha <= 1)
  m = m0 + m1
  if (m == 0L) return(list(p = numeric(0), N = NA_integer_, is_true = logical(0), U = numeric(0)))
  
  is_true = rep(FALSE, m)
  if (m0 > 0L) is_true[1:m0] = TRUE
  
  # Pre-generate U_1..U_m in S_j, plus U_{m+1} in (alpha,1]
  U = c(
    sapply(seq_len(m), function(j) runif(1L, (j - 1) * alpha / m, j * alpha / m)),
    runif(1L, alpha, 1)
  )
  
  # Draw N
  probs = .to_probs(m0, m1, alpha)
  N = sample.int(m + 1L, size = 1L, prob = probs)
  
  p = numeric(m)
  
  # False-null p-values: fixed 0's (independent of the true block)
  if (m1 > 0L) p[(m0 + 1L):m] = 0
  
  if (N <= m0) {
    # Case 1: 1 <= N <= m0 — put exactly N true nulls in S_N, rest > alpha
    if (m0 > 0L) {
      Tsel = sample.int(m0, N, replace = FALSE)
      p[Tsel] = U[N]
      rest = setdiff(seq_len(m0), Tsel)
      if (length(rest)) p[rest] = U[m + 1L]
    }
  } else if (N <= m) {
    # Case 2: m0+1 <= N <= m — put *all* true nulls in S_N
    if (m0 > 0L) p[1:m0] = U[N]
  } else {
    # Case 3: N = m+1 — put all true nulls > alpha
    if (m0 > 0L) p[1:m0] = U[m + 1L]
  }
  
  list(p = p, N = N, is_null = is_true, U = U)
}

# Guo Rao from n, pi0 and alpha only
gen_guo_rao_2 = function(n, pi0, alpha) {
  
  m = n
  m0 = floor(pi0 * m)
  m1 = m - m0
  
  sample_true_only(m0, m1, alpha)
  
}