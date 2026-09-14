library(posterior)
library(matrixStats)
library(splines)

# =============================================================================
# Ct latent-class model helpers
# =============================================================================


# =============================================================================
# Extract individual infection probabilities
# =============================================================================
# Returns a data frame with one row per fitted observation: median P(infection) and
# the observed ct, time since symptoms, and site.
get_prob_inf <- function(fit, stan_data, cts) {
  pc  <- exp(draws_of(as_draws_rvars(fit$draws("pC"))$pC))  # S x N x 2
  pii <- pc[, , 1] / (pc[, , 1] + pc[, , 2])
  data.frame(
    i    = seq_len(stan_data$N),
    pi   = colMedians(pii),
    ct   = cts$ct,
    t    = cts$Tsymp_hosp,
    site = cts$site
  )
}


# =============================================================================
# Posterior predictive density of Ct
# =============================================================================
# Returns a data frame with columns ct, med, ciL, ciU (posterior median and 95% CI
# of the predictive density, evaluated on a grid from `from` to `to`).
get_ppc_density <- function(fit, n_grid = 512, from = 10, to = 48, bw = 0.05) {
  y_rep_mat <- exp(as_draws_matrix(fit$draws("y_rep")))
  n_iter    <- nrow(y_rep_mat)
  den <- matrix(NA, ncol = n_grid, nrow = n_iter)
  for (i in seq_len(n_iter)) den[i, ] <- density(y_rep_mat[i, ], bw = bw, from = from, to = to)$y
  denq <- data.frame(ct = density(y_rep_mat[1, ], bw = bw, from = from, to = to)$x)
  denq[, c("med", "ciL", "ciU")] <- colQuantiles(den, probs = c(0.5, 0.025, 0.975))
  denq
}


# =============================================================================
# Population-average spline trajectory
# =============================================================================
# Computes the posterior mean Ct trajectory over a time grid, marginalised over
# the covariate distribution in the data. Returns a data frame with columns:
#   t, mu_mean, mu_lo, mu_hi  (posterior mean and 95% CI of the mean curve)
#   pred_lo, pred_hi           (80% posterior predictive interval)
# All Ct values are on the original (untransformed) scale.
get_spline_grid <- function(fit, stan_data, t_scaled, df_spline, t_grid = seq(0, 30, by = 0.5)) {
  draws_mat  <- as_draws_matrix(fit$draws(variables = c("c0", "beta", "alpha", "tau", "sdI")))
  c0_draws   <- draws_mat[, "c0"]
  beta_draws <- draws_mat[, grepl("^beta\\[",  colnames(draws_mat)), drop = FALSE]
  tau_draws  <- draws_mat[, grepl("^tau\\[",   colnames(draws_mat)), drop = FALSE]
  alpha_mat  <- draws_mat[, grepl("^alpha\\[", colnames(draws_mat)), drop = FALSE]

  S      <- nrow(draws_mat)
  NCov   <- stan_data$NCov
  NG     <- stan_data$NG
  maxNG  <- max(NG)
  covmat <- stan_data$cov

  alpha_arr <- array(0, dim = c(S, NCov, maxNG))
  for (c in seq_len(NCov)) {
    for (g in seq_len(NG[c])) {
      alpha_arr[, c, g] <- alpha_mat[, paste0("alpha[", c, ",", g, "]")]
    }
  }

  # Evaluate spline on a regular time grid (B is in data order, not time order).
  # Reuse the fitted basis via predict() so the knots match the ones `beta` was
  # estimated against; rebuilding with bs() takes its boundary knots from the
  # grid's own range, which only coincides with the data's by luck.
  t_grid_scaled <- (t_grid - attr(t_scaled, "scaled:center")) / attr(t_scaled, "scaled:scale")
  B_grid        <- if (inherits(stan_data$B, "bs")) {
    predict(stan_data$B, as.numeric(t_grid_scaled))
  } else {
    bs(as.numeric(t_grid_scaled), df = df_spline, degree = 4, intercept = TRUE)
  }

  mu_base_grid_draws <- sweep(beta_draws %*% t(B_grid), 1, c0_draws, "+")
  n_grid <- ncol(mu_base_grid_draws)

  # Population-average covariate offset (weighted by observed group frequencies)
  w <- matrix(0, nrow = NCov, ncol = maxNG)
  for (c in seq_len(NCov)) {
    for (g in seq_len(NG[c])) w[c, g] <- sum(covmat[, c] == g) / stan_data$N
  }
  delta_avg_draws <- vapply(seq_len(S), function(s) {
    sum(vapply(seq_len(NCov), function(c) tau_draws[s, c] * sum(w[c, seq_len(NG[c])] * alpha_arr[s, c, seq_len(NG[c])]), numeric(1)))
  }, numeric(1))

  mu_pop_grid_draws <- sweep(mu_base_grid_draws, 1, delta_avg_draws, "+")
  mu_pop_mean <- colMeans(mu_pop_grid_draws)
  mu_pop_ci   <- apply(mu_pop_grid_draws, 2, quantile, probs = c(0.025, 0.975))

  sigma_draws       <- draws_mat[, "sdI"]
  y_pred_grid_draws <- matrix(NA_real_, nrow = S, ncol = n_grid)
  for (s in seq_len(S)) {
    y_pred_grid_draws[s, ] <- rnorm(n_grid, mean = mu_pop_grid_draws[s, ], sd = sigma_draws[s])
  }
  pred_ci <- apply(y_pred_grid_draws, 2, quantile, c(0.1, 0.9))

  data.frame(
    t       = t_grid,
    mu_mean = exp(mu_pop_mean),
    mu_lo   = exp(mu_pop_ci[1, ]),
    mu_hi   = exp(mu_pop_ci[2, ]),
    pred_lo = exp(pred_ci[1, ]),
    pred_hi = exp(pred_ci[2, ])
  )
}


# =============================================================================
# Calculate P(infection | Ct<40) by covariate group
# =============================================================================
# Posterior P(true infection | Ct < 40) averaged within each covariate group,
# overall and by site. Groups with fewer than 2 observations below the threshold
# are left as NA. Returns a list of two data frames:
#   probI     - var, med, ciL, ciU, covar  (pooled across sites)
#   probIsite - var, site, med, ciL, ciU, covar
extract_pI_group <- function(fit, data, covlabs, sites, covars) {

  wb40      <- which(exp(data$y) < 40)
  n_sites   <- length(sites)
  NG        <- data$NG

  # Tabulate over the *full* set of levels, not just those surviving the Ct<40
  # subset: the g = 1..NG[c] loops below index these tables positionally, so a
  # level with no observations under the threshold would otherwise shift every
  # subsequent group onto the wrong label.
  NperG     <- lapply(seq_len(data$NCov), function(c)
    table(factor(data$cov[wb40, c], levels = seq_len(NG[c]))))
  NperGSite <- lapply(seq_len(data$NCov), function(c)
    table(factor(data$cov[wb40, c], levels = seq_len(NG[c])),
          factor(data$site[wb40],   levels = seq_len(n_sites))))

  pc   <- draws_of(as_draws_rvars(fit$draws("pC"))$pC)
  pc   <- exp(pc)
  pIi  <- pc[, , 1] / apply(pc, MARGIN = c(1, 2), sum)
  pIi  <- pIi[, wb40]
  n_iter    <- nrow(pIi)

  piG     <- vector("list", data$NCov)
  piGsite <- vector("list", data$NCov)

  for (covi in seq_len(data$NCov)) {
    tmp  <- array(NA, dim = c(n_iter, NG[covi]))
    tmpS <- array(NA, dim = c(n_iter, NG[covi], n_sites))
    for (g in seq_len(NG[covi])) {
      if (NperG[[covi]][g] > 1)
        tmp[, g] <- rowSums(pIi[, which(data$cov[wb40, covi] == g), drop = FALSE]) / NperG[[covi]][g]
      for (s in seq_len(n_sites)) {
        if (NperGSite[[covi]][g, s] > 1)
          tmpS[, g, s] <- rowSums(pIi[, which(data$cov[wb40, covi] == g & data$site[wb40] == s), drop = FALSE]) / NperGSite[[covi]][g, s]
      }
    }
    piG[[covi]]     <- tmp
    piGsite[[covi]] <- tmpS
  }

  probs     <- vector("list", data$NCov)
  probsSite <- vector("list", data$NCov)

  for (covi in seq_len(data$NCov)) {
    tmp <- data.frame(var = covlabs[[covi]], med = NA, ciL = NA, ciU = NA)
    for (g in seq_len(NG[covi])) {
      if (NperG[[covi]][g] > 1)
        tmp[g, 2:4] <- quantile(piG[[covi]][, g], c(0.5, 0.025, 0.975))
    }
    tmpS <- data.frame(var  = covlabs[[covi]],
                       site = sort(rep(sites, NG[covi])),
                       med  = NA, ciL = NA, ciU = NA)
    for (g in seq_len(NG[covi])) for (s in seq_len(n_sites)) {
      if (NperGSite[[covi]][g, s] > 1)
        tmpS[tmpS$var == covlabs[[covi]][g] & tmpS$site == sites[s], 3:5] <-
          quantile(piGsite[[covi]][, g, s], c(0.5, 0.025, 0.975))
    }
    probs[[covi]]           <- tmp
    probs[[covi]]$var       <- factor(probs[[covi]]$var, levels = covlabs[[covi]])
    probs[[covi]]$covar     <- covars[covi]
    probsSite[[covi]]       <- tmpS
    probsSite[[covi]]$covar <- covars[covi]
  }

  list(
    probI     = do.call("rbind", probs),
    probIsite = do.call("rbind", probsSite)
  )
}


# =============================================================================
# Extract covariate effects on Ct
# =============================================================================
# Posterior summaries of the scaled covariate effects on log Ct (Stan's
# alpha_eff = alpha * tau, i.e. effects on the log scale after centring within
# covariate). Returns a data frame with columns var, group, med, ciL, ciU,
# one row per covariate level, with `group` a factor in covlabs order.
get_covar_eff <- function(fit, data, covars, covlabs) {
  alpha <- draws_of(as_draws_rvars(fit$draws("alpha_eff"))$alpha_eff)
  L <- lapply(seq_along(covars), function(covi) {
    tmp <- data.frame(var = covars[covi], group = covlabs[[covi]], med = NA, ciL = NA, ciU = NA)
    for (g in seq_len(data$NG[covi]))
      tmp[g, 3:5] <- quantile(alpha[, covi, g], c(0.5, 0.025, 0.975))
    tmp$group <- factor(tmp$group, levels = covlabs[[covi]])
    tmp
  })
  do.call("rbind", L)
}


# =============================================================================
# Sensitivity and specificity at different Ct cut-offs
# =============================================================================
# cts_df: the data frame of Ct-positive records used to fit the model (rows aligned with data$y)
# sites:  sites to include in the ROC
get_roc <- function(fit, data, cutoffs, negnoct, cts_df, sites) {
  wroc    <- which(cts_df$site %in% sites)
  N       <- length(wroc)
  cts_val <- exp(data$y[wroc])

  pc  <- draws_of(as_draws_rvars(fit$draws("pC"))$pC)
  pc  <- exp(pc[, wroc, ])
  pIi <- pc[, , 1] / apply(pc, MARGIN = c(1, 2), sum)
  n_iter  <- nrow(pIi)

  true <- matrix(NA, nrow = n_iter, ncol = N)
  for (i in seq_len(N)) true[, i] <- rbinom(n_iter, size = 1, pIi[, i])

  roc <- lapply(cutoffs, function(co) {
    tmp <- data.frame(iter = seq_len(n_iter), tpr = NA, fpr = NA, tp = NA, fp = NA, fn = NA, cutoff = co)
    pos <- cts_val < co
    for (i in seq_len(n_iter)) {
      tp_i       <- sum( pos &  true[i, ] == 1)
      fp_i       <- sum( pos &  true[i, ] == 0)
      fn_i       <- sum(!pos &  true[i, ] == 1)
      tn_i       <- sum(!pos &  true[i, ] == 0) + negnoct
      tmp$tpr[i] <- tp_i / (tp_i + fn_i)
      tmp$fpr[i] <- fp_i / (fp_i + tn_i)
      tmp$tp[i]  <- tp_i
      tmp$fp[i]  <- fp_i
      tmp$fn[i]  <- fn_i
    }
    tmp
  })

  rocc <- do.call("rbind", roc)
  rocq <- data.frame(cutoff = cutoffs,
                     tpr = NA, tprL = NA, tprU = NA,
                     fpr = NA, fprL = NA, fprU = NA,
                     tp  = NA, tpL  = NA, tpU  = NA,
                     fp  = NA, fpL  = NA, fpU  = NA,
                     fn  = NA, fnL  = NA, fnU  = NA)
  for (i in seq_along(cutoffs)) {
    rocq[i,  2:4]  <- quantile(roc[[i]]$tpr, c(0.5, 0.025, 0.975))
    rocq[i,  5:7]  <- quantile(roc[[i]]$fpr, c(0.5, 0.025, 0.975))
    rocq[i,  8:10] <- quantile(roc[[i]]$tp,  c(0.5, 0.025, 0.975))
    rocq[i, 11:13] <- quantile(roc[[i]]$fp,  c(0.5, 0.025, 0.975))
    rocq[i, 14:16] <- quantile(roc[[i]]$fn,  c(0.5, 0.025, 0.975))
  }
  list(rocc = rocc, rocq = rocq)
}


# =============================================================================
# Population-marginalised infection / environmental density
# =============================================================================
# Averages the truncated-Normal density over all N individual mu[n]
# values within each posterior draw, giving the true population-marginal density.
get_marginal_dens <- function(fit, data) {
  mu_draws  <- as_draws_matrix(fit$draws("mu"))
  sdI_draws <- as_draws_matrix(fit$draws("sdI"))
  muE_draws <- as_draws_matrix(fit$draws("muE"))
  sdE_draws <- as_draws_matrix(fit$draws("sdE"))
  piI_draws <- as_draws_matrix(fit$draws("piI"))

  S          <- nrow(mu_draws)
  ct_grid    <- data$ct_range
  log_thresh <- log(if (!is.null(data$ct_threshold)) data$ct_threshold else 40)
  above_thresh <- ct_grid > log_thresh

  piI_avg <- as.numeric(piI_draws %*% data$NperSite) / sum(data$NperSite)

  inf_mat <- env_mat <- matrix(0, nrow = S, ncol = length(ct_grid))
  for (s in seq_len(S)) {
    mu_s  <- as.numeric(mu_draws[s, ])
    sdI_s <- as.numeric(sdI_draws[s, 1])
    muE_s <- as.numeric(muE_draws[s, 1])
    sdE_s <- as.numeric(sdE_draws[s, 1])
    dens_mat  <- outer(ct_grid, mu_s, function(ct, mu) dnorm(ct, mu, sdI_s))
    norms     <- pnorm(log_thresh, mu_s, sdI_s)
    dens_norm <- t(t(dens_mat) / norms)
    inf_dens  <- rowMeans(dens_norm)
    inf_dens[above_thresh] <- 0
    inf_mat[s, ] <- inf_dens
    env_mat[s, ] <- dnorm(ct_grid, muE_s, sdE_s)
  }

  pI_mat   <- sweep(inf_mat, 1, piI_avg,       "*")
  pEnv_mat <- sweep(env_mat, 1, 1.0 - piI_avg, "*")
  pAll_mat <- pI_mat + pEnv_mat

  summarise_dens <- function(mat) {
    data.frame(
      ct  = exp(ct_grid),
      med = colMedians(mat),
      ciL = colQuantiles(mat, probs = 0.025),
      ciU = colQuantiles(mat, probs = 0.975)
    )
  }
  list(pI = summarise_dens(pI_mat), pEnv = summarise_dens(pEnv_mat), pAll = summarise_dens(pAll_mat))
}


# =============================================================================
# Serology mixture model helpers
# =============================================================================


# =============================================================================
# Extract serology mixture model params
# =============================================================================
# Posterior median and 95% CrI for the serology mixture parameters.
# Returns a data frame with columns par, med, ciL, ciU and one row per
# parameter (sero, mu0, mu1, sd0, sd1).
extract_mixture_pars <- function(draws) {
  drawsdf <- as.data.frame(draws)
  pars <- c("sero", "mu0", "mu1", "sd0", "sd1")
  parests <- data.frame(par = pars, med = NA, ciL = NA, ciU = NA)
  for (p in seq_along(pars)) {
    parests[p, 2:4] <- quantile(drawsdf[, pars[p]], c(0.5, 0.025, 0.975))
  }
  parests
}

# =============================================================================
# Extract serology mixture model fits
# =============================================================================
# Fitted negative and positive component densities evaluated on sdata$y_fit.
# Returns a long data frame with columns component ("Negative"/"Positive"),
# log_titer, med, ciL, ciU.
extract_mixture_fits <- function(fit, sdata) {
  dists  <- c("fitNeg", "fitPos")
  labels <- c("Negative", "Positive")
  fits <- lapply(seq_along(dists), function(i) {
    mat <- fit$draws(dists[i], format = "matrix")
    tmp <- data.frame(component = labels[i], log_titer = sdata$y_fit, med = NA, ciL = NA, ciU = NA)
    tmp[, 3:5] <- t(apply(mat, 2, quantile, probs = c(0.5, 0.025, 0.975)))
    tmp
  })
  fits <- do.call(rbind, fits)
  rownames(fits) <- NULL
  fits
}


# =============================================================================
# Plot serology mixture model fits
# =============================================================================
# Histogram of observed log(MFI) overlaid with the fitted mixture components
# and their 95% CrI ribbons. Returns a ggplot object; `antigen` is used as the
# panel subtitle.
plot_mixture_fits <- function(fit, sdata, antigen) {
  titer_fit <- extract_mixture_fits(fit, sdata)
  ggplot() +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab("log(MFI)") +
    geom_histogram(
      data = data.frame(titer = sdata$y),
      aes(x = titer, y = after_stat(density)),
      fill = "grey70", bins = 25, color = "white"
    ) +
    geom_line(data = titer_fit, aes(x = log_titer, y = med, color = component), linewidth = 1) +
    geom_ribbon(data = titer_fit, aes(x = log_titer, ymin = ciL, ymax = ciU, fill = component), alpha = 0.2) +
    labs(subtitle = antigen)
}

# =============================================================================
# Extract mixture model serostatus results
# =============================================================================
# Per-sample posterior probability of belonging to the seropositive component,
# computed as mean(exp(pC[2,n] - log_lik[n])) across draws. Returns a data frame
# with columns id, visit, log_mfi, prob_pos and seropos (prob_pos >= threshold).
extract_serostatus <- function(fit, dfa, threshold = 0.5) {
  draws <- fit$draws(c("pC", "log_lik"), format = "matrix")
  N     <- nrow(dfa)

  prob_pos <- vapply(seq_len(N), function(n) {
    mean(exp(draws[, paste0("pC[2,", n, "]")] - draws[, paste0("log_lik[", n, "]")]))
  }, numeric(1))

  data.frame(
    id       = dfa$id,
    visit    = dfa$visit,
    log_mfi  = log(dfa$mfi),
    prob_pos = prob_pos,
    seropos  = prob_pos >= threshold
  )
}

# =============================================================================
# Generate chain initialization values for serology mixture model
# =============================================================================
# Random starting values for the serology mixture model, one list per chain.
# mu1 is initialised above mu0 to discourage label switching between the
# seronegative and seropositive components.
make_sero_inits <- function(n_chains) {
  lapply(seq_len(n_chains), function(i) {
    list(
      sero = runif(1, 0.3, 0.7),
      mu0  = runif(1, 2.5, 5),
      mu1  = runif(1, 4, 6),
      sd0  = runif(1, 0.2, 1),
      sd1  = runif(1, 0.2, 1)
    )
  })
}
