# =============================================================================
# Fit latent class model to MPXV Ct value data, DRC
# =============================================================================
# Two-component mixture (true infection vs environmental contamination) over
# log Ct values, with a spline in time since symptom onset and covariate
# effects. Writes all model outputs to results/main; MakeFigures.R reads them.

library(cmdstanr)
library(splines)
library(ggridges)
library(bayesplot)
library(posterior)
library(GGally)
library(loo)
library(stringr)
library(matrixStats)
library(cowplot)
library(ggh4x)

# =============================================================================
# Setup
# =============================================================================

# ===== Paths =====
# Run this script from the repo root, or set repo_dir explicitly.
repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")
stan_dir <- file.path(repo_dir, "stan")
out_dir  <- file.path(repo_dir, "results/main")
dir.create(out_dir, showWarnings = FALSE)


# ===== Helper functions =====
source(file.path(repo_dir, "Utils.R"))


source(file.path(repo_dir, "LoadCaseCts.R"))

# ===== Data =====
df  <- read_case_cts(data_dir)
env <- read.csv(file.path(data_dir, "EnvironmentalSamplingCts.csv"))

noamp <- df[df$ct == 0, ]   # qPCR results with no amplification
cts   <- df[df$ct > 0, ]    # qPCR results with Ct values

sites <- levels(factor(cts$site))


# ===== Covariates =====
covars <- c("sex", "sample_type", "age_group", "sex_contact",
            "recent_vacc", "hiv", "lesions_anogen", "any_lymph", "lesions_palmsole")

# covar matrix & labels
# Codes and labels come from the SAME factor call so that code g always maps to
# covlabs[[v]][g]. Going via as.matrix() would coerce every column to character
# (which cmdstanr rejects) and re-sort factor levels alphabetically, silently
# mislabelling covariates whose level order is not alphabetical (e.g. age_group).
# read_case_cts() is what puts age_group back in its natural order after the CSV
# read; see LoadCaseCts.R.
NCov    <- length(covars)
covlabs <- lapply(covars, function(v) levels(factor(cts[[v]])))
names(covlabs) <- covars
NG      <- vapply(covlabs, length, integer(1))
cov     <- vapply(covars,
                  function(v) as.integer(factor(cts[[v]], levels = covlabs[[v]])),
                  integer(nrow(cts)))

# indices for model input
ind_cov      <- matrix(1, nrow = NCov, ncol = 2)
ind_cov[, 2] <- cumsum(NG)
ind_cov[, 1] <- cumsum(NG) - NG + 1


# ===== Spline basis =====
t_scaled  <- scale(cts$Tsymp_hosp, center = TRUE, scale = TRUE)
df_spline <- 5
B <- bs(as.numeric(t_scaled), df = df_spline, degree = 4, intercept = TRUE)
K <- ncol(B)


# ===== Ct threshold =====
ct_thresh <- 40


# ===== Stan data =====
stan_data <- list(
  N            = nrow(cts),
  NE           = nrow(env),
  N_sites      = length(sites),
  y            = log(cts$ct),
  t            = cts$Tsymp_hosp,
  yE           = log(env$envct),
  ct_range     = log(seq(1, 45)),
  site         = as.numeric(factor(cts$site)),
  NperSite     = as.numeric(table(cts$site)),
  N_ct_thresh  = as.numeric(table(factor(cts$site[cts$ct < ct_thresh], levels = sites))),
  ct_threshold = ct_thresh,
  NCov     = NCov,
  NG       = NG,
  cov      = cov,
  ind_cov  = ind_cov,
  K        = K,
  B        = B
)

# ===== Chain initialisation =====
init_chains <- function(stan_data, nChains) {
  lapply(seq_len(nChains), function(i) {
    list(
      piI        = runif(stan_data$N_sites, 0.7, 0.9),
      muE        = runif(1, 3.6, 3.9),
      sdE        = runif(1, 0.02, 0.1),
      sdI        = runif(1, 0.05, 0.2),
      c0         = runif(1, 3.3, 3.7),
      beta       = runif(stan_data$K, -2, 2),
      sigma_beta = runif(1, 1, 2),
      alpha_raw  = runif(sum(stan_data$NG), -0.5, 0.5),
      tau        = runif(stan_data$NCov, 0, 1)
    )
  })
}
# Fixed seed so the fit and every stochastic quantity derived from it are
# reproducible. Four separate sources of randomness have to be pinned, not
# just the sampler: the chain inits below (runif), Stan's own sampler, the
# rnorm() predictive draws in get_spline_grid(), and the rbinom() true-status
# draws in get_roc() that produce the ROC curve and the case counts. Each
# stochastic step is seeded explicitly rather than relying on one set.seed()
# at the top, so that reordering the script cannot silently change results.
main_seed <- 20260821

set.seed(main_seed)                     # governs init_chains() below
ini <- init_chains(stan_data, nChains = 3)


# =============================================================================
# Compile and fit
# =============================================================================
check_cmdstan_toolchain()
# set_cmdstan_path("path/to/cmdstan-X.X.X")
mod <- cmdstan_model(file.path(stan_dir, "LatentClassCtModel.stan"), pedantic = FALSE)

# fit model
fit <- mod$sample(
  data            = stan_data,
  chains          = 3,
  parallel_chains = 3,
  iter_sampling   = 2000,
  iter_warmup     = 1000,
  refresh         = 10,
  seed            = main_seed,
  init            = ini
)


# =============================================================================
# Convergence checks
# =============================================================================
# Pull only the scalar / low-dimensional parameters. fit$draws(format = "df")
# without `variables` also materialises mu, pC, log_lik and y_rep
# (~13,800 columns x 6,000 draws, ~0.7 GB) purely to draw a few traceplots.
chains <- fit$draws(
  variables = c("lp__", "piI", "truepI", "truepI_all", "muE", "sdE", "sdI",
                "c0", "beta", "sigma_beta", "alpha_raw", "tau"),
  format    = "df"
)
color_scheme_set("mix-blue-red")

# traceplots
mcmc_trace(chains, regex_pars = c("piI", "truepI", "lp__"))
mcmc_trace(chains, regex_pars = c("alpha_raw", "tau"))
mcmc_trace(chains, pars = c("muE", "sdE", "sdI"))
mcmc_trace(chains, regex_pars = c("c0", "beta"))

# rhats
rhat_check <- fit$summary(
  variables = c("piI", "muE", "sdE", "sdI", "c0", "sigma_beta", "tau", "alpha_raw"),
  "rhat"
)
print(rhat_check[, c("variable", "rhat")])


# =============================================================================
# Posterior summaries
# =============================================================================

# ===== Extract draws =====
cc <- as.data.frame(chains)


# ===== LOO / WAIC =====
ll      <- fit$draws("log_lik")
loo_1   <- loo(ll)
waicloo <- rbind(as.data.frame(loo_1$estimates), as.data.frame(waic(ll)$estimates))


# ===== Parameter summaries =====
cc_cols <- c(paste0("truepI[", seq_len(stan_data$N_sites), "]"), "truepI_all", "muE", "sdE", "sdI")
param_summary <- data.frame(
  pars = c(rep("probI", stan_data$N_sites), "probI_all", "muE", "sdE", "sdI"),
  site = c(sites, NA, NA, NA, NA),
  t(sapply(cc_cols, function(p) quantile(cc[[p]], c(med = 0.5, ciL = 0.025, ciU = 0.975))))
)
param_summary


# =============================================================================
# Derived quantities
# =============================================================================

# ===== Individual infection probabilities =====
pii <- get_prob_inf(fit, stan_data, cts)
ggplot(pii, aes(ct, pi)) + geom_point()


# ===== Posterior predictive check =====
# Overall Ct distribution.
denq <- get_ppc_density(fit)

fitplot <- ggplot() +
  geom_histogram(data = cts, aes(ct, after_stat(density)), bins = 30, col = "grey40", fill = "grey") +
  theme_bw() +
  geom_line(data = denq, aes(ct, med), col = "blue") +
  geom_ribbon(data = denq, aes(ct, ymin = ciL, ymax = ciU), alpha = 0.5, fill = "dodgerblue2") +
  theme(text = element_text(size = 14)) + xlab("Ct value")
fitplot


# ===== Spline trajectory =====
# Population-average Ct vs time since symptoms.
set.seed(main_seed + 1)                 # get_spline_grid() draws predictive rnorm()
df_grid <- get_spline_grid(fit, stan_data, t_scaled, df_spline)


# ===== Infection / environmental densities =====
# True population marginal.
dens <- get_marginal_dens(fit, stan_data)
pdfs <- do.call(rbind, Map(function(df, g) { df$group <- g; df }, dens, names(dens)))


# ===== Covariate effects on Ct =====
# Pretty-print the French HIV codes. Mapped BY NAME, not by position: hiv is the
# only covariate with mixed-case levels, so its ordering is collation-dependent
# (LC_COLLATE=C gives Negatif < Positif < missing; en_US.UTF-8 gives missing
# first). A positional override is silently wrong under the other locale.
hiv_labels <- c(Negatif = "negative", Positif = "positive", missing = "missing")
stopifnot(all(covlabs[["hiv"]] %in% names(hiv_labels)))
covlabs[["hiv"]] <- unname(hiv_labels[covlabs[["hiv"]]])
alpha <- get_covar_eff(fit, stan_data, covars, covlabs)

alpha$var <- factor(alpha$var,
  levels = c("age_group", "sex", "sample_type", "sex_contact", "hiv", "recent_vacc", "any_lymph", "lesions_anogen", "lesions_palmsole"),
  labels = c("Age", "Sex", "Sample type", "Sexual contact", "HIV", "Vaccinated", "Lymphadenopathy", "Genital/perianal lesions", "Palm/sole lesions")
)


# ===== Infection probability by covariate group =====
picov <- extract_pI_group(fit, stan_data, covlabs = covlabs,
                          sites = sites, covars = covars)
probI <- picov$probI
probI <- probI[!is.na(probI$med), ]

probI$covar <- factor(probI$covar,
  levels = c("age_group", "sex", "sample_type", "sex_contact", "hiv", "recent_vacc", "any_lymph", "lesions_anogen", "lesions_palmsole"),
  labels = c("Age", "Sex", "Sample type", "Sexual contact", "HIV", "Vaccinated", "Lymphadenopathy", "Genital/perianal lesions", "Palm/sole lesions")
)


# =============================================================================
# ROC analysis
# =============================================================================
set.seed(main_seed + 2)                 # get_roc() draws true status via rbinom()
ROC <- get_roc(fit, stan_data, cutoffs = seq(20, ct_thresh),
               negnoct = nrow(noamp[noamp$site %in% sites, ]), cts_df = cts, sites = sites)

# quick plot
rocplot <- ggplot(ROC$rocq[ROC$rocq$cutoff > 29, ], aes(fpr, tpr, label = cutoff)) +
  theme_bw() +
  geom_line(data = ROC$rocc, aes(fpr, tpr, group = iter), alpha = 0.1, col = "grey80") +
  xlim(0, NA) +
  geom_linerange(aes(ymin = tprL, ymax = tprU), col = "blue") +
  geom_point(col = "blue") +
  geom_linerange(aes(xmin = fprL, xmax = fprU), col = "blue") +
  geom_text(vjust = 1.5, hjust = 1.3, col = "indianred2") +
  theme(text = element_text(size = 14), legend.position = "none") +
  xlab("False positive rate") + ylab("Sensitivity")
rocplot

# calculate youdens J statistic
ROC$rocq$youden <- ROC$rocq$tpr + (1 - ROC$rocq$fpr) - 1


# =============================================================================
# Save results
# =============================================================================
fit$save_object(file.path(out_dir, "Fit.RDS"))  # saveRDS() leaves the fit pointing at a tempdir that is gone next session
saveRDS(stan_data, file.path(out_dir, "Data.RDS"))
saveRDS(ROC,       file.path(out_dir, "ROC.RDS"))
write.csv(waicloo,                            file.path(out_dir, "WAICLOO.csv"))
write.csv(df_grid,                            file.path(out_dir, "TimeSinceSymps.csv"), row.names = FALSE)
write.csv(pdfs,                               file.path(out_dir, "FitPDFs.csv"),        row.names = FALSE)
write.csv(param_summary,                      file.path(out_dir, "Params.csv"),         row.names = FALSE)
write.csv(probI,                              file.path(out_dir, "pI_covars.csv"),      row.names = FALSE)
write.csv(alpha,                              file.path(out_dir, "covar_coef.csv"),     row.names = FALSE)
write.csv(pii[, c("pi", "ct", "t", "site")],  file.path(out_dir, "pii.csv"),            row.names = FALSE)

