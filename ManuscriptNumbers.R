# =============================================================================
# Manuscript numbers audit
# =============================================================================
# Recomputes, from the public data and saved model outputs in this repository,
# every number reported in O'Driscoll et al. (Lancet Infect Dis, 2026),
# "The risk of mpox false positive results in high transmission settings",
# and checks each one against the value printed in the accepted proofs
# (26TLID0924_Azman, fourth proof, saved 17:55 09-Sept-2026).
#
# Usage (from the repo root):
#     source("ManuscriptNumbers.R")
#   or
#     Rscript ManuscriptNumbers.R
#
# Every line printed carries one of three verdicts:
#   MATCH      recomputed value agrees with the proofs
#   DIFFERS    recomputed value disagrees -- inspect before the proofs go back
#   NO DATA    input needed is not in the public repo (see the `note` column)
#
# A section near the end ("EDITOR QUERIES") answers the [A: ...] author queries
# in this proof that ask for a number.
#
# A machine-readable copy of the audit is written to
# results/ManuscriptNumbers.csv.
#
# Base R only; no model refitting. Credible intervals for the relative risks in
# figure 2 need the posterior draws (results/<dir>/Fit.RDS, not committed --
# run FitMainModel.R). Without it the script still prints RR point estimates.
# =============================================================================


# =============================================================================
# Configuration
# =============================================================================

repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")

# Directory holding the fitted-model outputs written by FitMainModel.R.
# Override with the environment variable MPOX_RESULTS_DIR.
results_dir <- Sys.getenv("MPOX_RESULTS_DIR", unset = "results/main")

sero_file <- file.path(repo_dir, "results/serology/SeroResults.csv")
audit_out <- file.path(repo_dir, "results/ManuscriptNumbers.csv")

# Recompute relative-risk credible intervals from the posterior draws when
# Fit.RDS is available. Costs a few minutes and several GB of RAM; the result is
# cached as RR_covars.csv in results_dir and reused on later runs.
use_fit_draws <- TRUE


# =============================================================================
# Audit bookkeeping
# =============================================================================

.audit <- data.frame(section = character(), item = character(),
                     published = character(), computed = character(),
                     verdict = character(), note = character(),
                     stringsAsFactors = FALSE)

W_ITEM <- 44; W_PUB <- 24; W_COMP <- 24

.section_now <- ""

sec <- function(title) {
  .section_now <<- title
  cat("\n", strrep("=", 122), "\n", sep = "")
  cat(title, "\n")
  cat(strrep("=", 122), "\n")
  cat(sprintf("  %-*s  %-*s  %-*s  %s\n",
              W_ITEM, "", W_PUB, "PUBLISHED", W_COMP, "RECOMPUTED", "VERDICT"))
}

add <- function(item, published, computed, verdict, note = "") {
  .audit <<- rbind(.audit, data.frame(
    section = .section_now, item = item, published = published,
    computed = computed, verdict = verdict, note = note,
    stringsAsFactors = FALSE))
  cat(sprintf("  %-*s  %-*s  %-*s  %s%s\n",
              W_ITEM, substr(item, 1, W_ITEM),
              W_PUB,  substr(published, 1, W_PUB),
              W_COMP, substr(computed, 1, W_COMP),
              verdict,
              if (nzchar(note)) paste0("  <- ", note) else ""))
  invisible(NULL)
}

# Auto-verdict numeric check. `pub` and `comp` are numeric vectors of the same
# length; they agree when every element is within `tol` (recycled). `style` sets
# the rendering:
#   "ci"    point estimate then interval  -> 35.0% (31.0-39.0)
#   "range" a range                       -> 17-27
#   "npct"  a count and its percentage    -> 1278 (46.9%)
num <- function(item, pub, comp, tol = 0.5, digits = 1, unit = "", note = "",
                style = c("ci", "range", "npct")) {
  style <- match.arg(style)
  f <- function(x) {
    if (all(is.na(x))) return("--")
    s <- formatC(x, format = "f", digits = digits)
    if (length(x) == 1) return(paste0(s, unit))
    switch(style,
      ci    = paste0(s[1], unit, " (", paste(s[-1], collapse = "-"), ")"),
      range = paste(s, collapse = "-"),
      npct  = sprintf("%d (%s%%)", round(x[1]), formatC(x[2], format = "f", digits = 1)))
  }
  keep <- !is.na(pub) & !is.na(comp)
  ok <- length(pub) == length(comp) && any(keep) &&
        all(abs(pub[keep] - comp[keep]) <= rep_len(tol, length(pub))[keep] + 1e-9)
  if (!all(keep))
    note <- trimws(paste(note, "(no CrI available; point estimate compared)"))
  add(item, f(pub), f(comp), if (isTRUE(ok)) "MATCH" else "DIFFERS", note)
}

# Count with its within-column percentage.
cnt <- function(item, pub_n, pub_pct, k, n, tol = c(0, 0.05), note = "")
  num(item, c(pub_n, pub_pct), c(k, pct(k, n)), tol = tol, style = "npct", note = note)

# Manual verdict, for values that are not a single number.
txt <- function(item, published, computed, verdict, note = "")
  add(item, published, computed, verdict, note)

# Number reported in the paper that the public data cannot reproduce.
gap <- function(item, published, why)
  add(item, published, "--", "NO DATA", why)

pct  <- function(k, n) k / n * 100
np   <- function(k, n, digits = 1)
  sprintf("%d (%s%%)", k, formatC(pct(k, n), format = "f", digits = digits))


# =============================================================================
# Load inputs
# =============================================================================

source(file.path(repo_dir, "LoadCaseCts.R"))

stopifnot(dir.exists(data_dir))
if (!dir.exists(results_dir))
  stop("results_dir not found: ", results_dir,
       "\nRun FitMainModel.R first, or set MPOX_RESULTS_DIR.")

df    <- read_case_cts(data_dir)
env   <- read.csv(file.path(data_dir, "EnvironmentalSamplingCts.csv"))
names(env)[1] <- "site"                 # strip the UTF-8 BOM on the header

cts   <- df[df$ct > 0, ]                # detectable OPXV DNA
noamp <- df[df$ct == 0, ]               # no amplification

sites   <- c("Goma", "Kamituga", "Kinshasa", "Uvira")
gx_sites <- c("Kamituga", "Uvira")      # GeneXpert; the other two used RADI

params    <- read.csv(file.path(results_dir, "Params.csv"))
pI_covars <- read.csv(file.path(results_dir, "pI_covars.csv"))
tss       <- read.csv(file.path(results_dir, "TimeSinceSymps.csv"))
waicloo   <- read.csv(file.path(results_dir, "WAICLOO.csv"), row.names = 1)
roc       <- readRDS(file.path(results_dir, "ROC.RDS"))$rocq
stan_data <- readRDS(file.path(results_dir, "Data.RDS"))

sero <- if (file.exists(sero_file)) read.csv(sero_file) else NULL

# One row per serology participant, with serostatus and the Ct-model class.
# Built here because both the abstract and the figure 3B section check them.
if (!is.null(sero)) {
  ss      <- sero[!duplicated(sero$id), ]
  seropos <- ss$statusV2 == "Positive"
  sero_fp <- ss$ct > 0 & ss$ct < 40 & ss$prob_inf <  0.5   # inferred false positives
  sero_tp <- ss$ct > 0 & ss$ct < 40 & ss$prob_inf >= 0.5   # inferred true positives
}

par_med <- function(p, s = NA) {
  r <- if (is.na(s)) params[params$pars == p, ] else params[params$pars == p & params$site %in% s, ]
  unlist(r[1, c("X50.", "X2.5.", "X97.5.")])
}
roc_at <- function(cut, cols) unlist(roc[roc$cutoff == cut, cols])


# =============================================================================
# Covariate labels: recover them from the saved Stan data
# =============================================================================
# extract_pI_group() writes one row per covariate group in group-index order,
# labelled with whatever level order the fitting session happened to produce.
# Runs made before the covariate-coding fix (August 2026) coded age_group
# and hiv in alphabetical order ("0-4" < "10-14" < ... < "5-9") while labelling
# them in natural order, which permutes the age labels in pI_covars.csv. The
# true label of group g is recoverable by cross-tabulating stan_data$cov against
# the raw data, so do that rather than trusting the CSV.

covars <- c("sex", "sample_type", "age_group", "sex_contact",
            "recent_vacc", "hiv", "lesions_anogen", "any_lymph", "lesions_palmsole")
pretty_covar <- c(sex = "Sex", sample_type = "Sample type", age_group = "Age",
                  sex_contact = "Sexual contact", recent_vacc = "Vaccinated",
                  hiv = "HIV", lesions_anogen = "Genital/perianal lesions",
                  any_lymph = "Lymphadenopathy", lesions_palmsole = "Palm/sole lesions")
hiv_labels <- c(Negatif = "negative", Positif = "positive", missing = "missing")

stopifnot(nrow(stan_data$cov) == nrow(cts))
pI_covars$covar_raw <- covars[match(pI_covars$covar, pretty_covar[covars])]

recovered <- lapply(seq_along(covars), function(i) {
  tb <- table(stan_data$cov[, i], cts[[covars[i]]])
  lab <- colnames(tb)[apply(tb, 1, which.max)]
  stopifnot(!anyDuplicated(lab))        # coding must be a bijection
  if (covars[i] == "hiv") lab <- unname(hiv_labels[lab])
  lab
})
names(recovered) <- covars

natural <- lapply(covars, function(v) {
  lab <- levels(factor(cts[[v]]))
  if (v == "hiv") lab <- unname(hiv_labels[lab]) else lab
})
names(natural) <- covars

# Flag covariates whose labels in pI_covars.csv are attached to the wrong group.
csv_labels <- lapply(covars, function(v) pI_covars$var[pI_covars$covar == pretty_covar[[v]]])
names(csv_labels) <- covars
relabelled <- vapply(covars, function(v)
  length(csv_labels[[v]]) == length(recovered[[v]]) &&
    !identical(as.character(csv_labels[[v]]), recovered[[v]]), logical(1))

# Print groups in the data's own level order, whatever order the model coded them in.
display_order <- lapply(covars, function(v) natural[[v]][natural[[v]] %in% recovered[[v]]])
names(display_order) <- covars

# Re-attach the recovered labels to pI_covars, block by block.
pI_covars$label <- NA_character_
for (v in covars) {
  w <- which(pI_covars$covar_raw == v)
  if (length(w) == length(recovered[[v]])) {
    pI_covars$label[w] <- recovered[[v]]
  } else {
    pI_covars$label[w] <- pI_covars$var[w]      # a group was dropped; keep as-is
    warning("pI_covars.csv has ", length(w), " rows for '", v,
            "' but the model had ", length(recovered[[v]]),
            " groups; labels for this covariate left unchanged.")
  }
}

pI <- function(v, g) {
  r <- pI_covars[pI_covars$covar_raw == v & pI_covars$label == g, ]
  if (nrow(r) != 1) return(c(NA, NA, NA))
  unlist(r[1, c("med", "ciL", "ciU")])
}


# =============================================================================
# Relative risks of a false positive result (figure 2)
# =============================================================================
# RR_g = P(false positive | Ct<40, group g) / P(false positive | Ct<40, ref).
# The credible interval is the interval of the ratio, so it has to be formed
# draw by draw; the marginal intervals in pI_covars.csv cannot be combined.

rr_ref <- c(sex = "female", sample_type = "lesion", age_group = "20-29",
            sex_contact = "no", recent_vacc = "no", hiv = "negative",
            lesions_anogen = "no", any_lymph = "no", lesions_palmsole = "no")

rr_cache <- file.path(results_dir, "RR_covars.csv")
fit_file <- file.path(results_dir, "Fit.RDS")

compute_rr_draws <- function() {
  if (!requireNamespace("posterior", quietly = TRUE) ||
      !requireNamespace("cmdstanr", quietly = TRUE)) {
    message("posterior and cmdstanr are needed to read Fit.RDS; skipping RR intervals.")
    return(NULL)
  }
  loadNamespace("cmdstanr")   # the fit is an R6 object of a cmdstanr class
  message("Reading posterior draws from ", fit_file,
          " (slow; cached to RR_covars.csv afterwards) ...")
  fit <- readRDS(fit_file)
  pc  <- exp(posterior::draws_of(posterior::as_draws_rvars(fit$draws("pC"))$pC))
  pIi <- pc[, , 1] / apply(pc, MARGIN = c(1, 2), sum)
  rm(pc); gc()

  wb40 <- which(exp(stan_data$y) < 40)
  pIi  <- pIi[, wb40, drop = FALSE]
  covm <- stan_data$cov[wb40, , drop = FALSE]

  out <- list()
  for (i in seq_along(covars)) {
    v    <- covars[i]
    labs <- recovered[[v]]
    grp  <- lapply(seq_along(labs), function(g) {
      w <- which(covm[, i] == g)
      if (length(w) < 2) NULL else rowMeans(pIi[, w, drop = FALSE])
    })
    names(grp) <- labs
    ref <- grp[[rr_ref[[v]]]]
    if (is.null(ref)) next
    for (g in labs) {
      if (is.null(grp[[g]])) next
      q <- quantile((1 - grp[[g]]) / (1 - ref), c(0.5, 0.025, 0.975))
      out[[length(out) + 1]] <- data.frame(
        covar = v, label = g, ref = rr_ref[[v]],
        n = sum(covm[, i] == which(labs == g)),
        rr = q[[1]], rrL = q[[2]], rrU = q[[3]], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# A cached RR table is only valid for the reference groups it was built with,
# so discard it when rr_ref has changed since (or when it predates the `ref`
# column). Silently reusing it would report ratios against the wrong baseline.
rr_tab <- NULL
if (file.exists(rr_cache)) {
  cached <- read.csv(rr_cache, stringsAsFactors = FALSE)
  same_ref <- "ref" %in% names(cached) &&
    identical(sort(unique(paste(cached$covar, cached$ref))),
              sort(unique(paste(names(rr_ref), unname(rr_ref)))))
  # A cache older than the fit was built from a previous run of the model.
  fresh <- !file.exists(fit_file) ||
    file.mtime(rr_cache) >= file.mtime(fit_file)
  if (same_ref && fresh) rr_tab <- cached
  else if (!same_ref) message("RR cache was built with different reference groups; recomputing.")
  else message("RR cache predates ", fit_file, "; recomputing.")
}
if (is.null(rr_tab) && use_fit_draws && file.exists(fit_file)) {
  rr_tab <- tryCatch(compute_rr_draws(),
                     error = function(e) { message("RR draws failed: ", conditionMessage(e)); NULL })
  if (!is.null(rr_tab)) write.csv(rr_tab, rr_cache, row.names = FALSE)
}

# Fall back to point estimates from the group medians when there are no draws.
rr_have_ci <- !is.null(rr_tab)
if (!rr_have_ci) {
  rr_tab <- do.call(rbind, lapply(covars, function(v) {
    ref <- pI(v, rr_ref[[v]])[1]
    do.call(rbind, lapply(recovered[[v]], function(g) {
      p <- pI(v, g)[1]
      data.frame(covar = v, label = g, ref = rr_ref[[v]],
                 n = sum(stan_data$cov[, match(v, covars)] ==
                           match(g, recovered[[v]])),
                 rr = (1 - p) / (1 - ref), rrL = NA_real_, rrU = NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
}

rr_of <- function(v, g) {
  r <- rr_tab[rr_tab$covar == v & rr_tab$label == g, ]
  if (nrow(r) != 1) return(c(NA, NA, NA))
  unlist(r[1, c("rr", "rrL", "rrU")])
}
rr_str <- function(x) {
  if (is.na(x[1])) return("--")
  if (is.na(x[2])) sprintf("%.2f (no CrI)", x[1])
  else sprintf("%.2f (%.2f-%.2f)", x[1], x[2], x[3])
}


# =============================================================================
cat("\n")
cat(strrep("#", 122), "\n")
cat("# MANUSCRIPT NUMBERS -- O'Driscoll et al., Lancet Infect Dis 2026 (proof 26TLID0924)\n")
cat("# repo:        ", repo_dir, "\n", sep = "")
cat("# model output:", results_dir, "\n", sep = "")
cat("# run:         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat(strrep("#", 122), "\n")

if (any(relabelled)) {
  cat("\n!! WARNING: group labels in ", results_dir, "/pI_covars.csv are attached to the wrong\n",
      "   groups for: ", paste(covars[relabelled], collapse = ", "), "\n", sep = "")
  cat("   That output predates the covariate-coding fix. This script uses the labels\n")
  cat("   recovered from Data.RDS, so its estimates are correctly labelled, but any\n")
  cat("   figure or text built straight from that CSV is not. Recovered order:\n")
  for (v in covars[relabelled])
    cat("     ", v, ": ", paste(recovered[[v]], collapse = " | "), "\n", sep = "")
}
if (!rr_have_ci)
  cat("\n!! NOTE: ", fit_file, " not found -- relative risks are point estimates with no CrI.\n", sep = "")


# =============================================================================
sec("SUMMARY / ABSTRACT")

num("Suspected cases with Ct value data", 2724, nrow(cts), tol = 0, digits = 0)
num("False positives, Ct<40 (%, 95% CrI)",
    c(35, 31, 39), 100 * (1 - par_med("probI_all")[c(1, 3, 2)]), tol = 0.5, unit = "%")
num("False positive rate, Ct<37 cutoff",
    c(18, 16, 20), 100 * roc_at(37, c("fpr", "fprL", "fprU")), tol = 0.5, unit = "%")
num("False positive rate, Ct<34 cutoff",
    c(5, 3, 6), 100 * roc_at(34, c("fpr", "fprL", "fprU")), tol = 0.5, unit = "%")

# Added to the abstract in this proof, flagged "[A: this needs to be added to
# the main Results]". The two denominators are the inferred false and true
# positives among the serology subset with a Ct value below 40.
if (!is.null(sero)) {
  cnt("Inferred FPs with no serological evidence", 31, 88.6,
      sum(sero_fp & !seropos), sum(sero_fp), tol = c(0, 0.5))
  cnt("Inferred true positives with serological evidence", 20, 66.7,
      sum(sero_tp & seropos), sum(sero_tp), tol = c(0, 0.5))
}


# =============================================================================
sec("METHODS")

gap("Environmental samples tested, Uvira", "41",
    "only the 11 OPXV-positive Uvira swabs are in EnvironmentalSamplingCts.csv")
gap("Environmental samples tested, Kamituga", "40",
    "only the 23 OPXV-positive Kamituga swabs are in EnvironmentalSamplingCts.csv")
if (!is.null(sero))
  num("Participants in serological follow-up", 118, length(unique(sero$id)), tol = 0, digits = 0)
gap("Cases in masked photo review", "116",
    "clinician review data (ratings, Likert scores) are not in this repo")
num("Age 20-29 reference group, n with Ct<40", 695,
    sum(cts$age_group == "20-29" & cts$ct < 40), tol = 0, digits = 0)
num("ROC Ct cutoffs evaluated (from-to)", c(20, 40), range(roc$cutoff),
    tol = 0, digits = 0, style = "range")
gap("Reporting period, per site and overall",
    "May 2, 2024-April 27, 2026",
    "CaseCts.csv carries no collection dates; the table 1 periods are not recomputable")
cat("  (Sampler settings -- 3 chains, 2000 iterations after 1000 warm-up -- are set in\n",
    "   FitMainModel.R and are not recoverable from the saved summaries.)\n", sep = "")


# =============================================================================
sec("RESULTS -- STUDY POPULATION (paragraph 1)")

num("Individuals with a qPCR result", 4004, nrow(df), tol = 0, digits = 0)
num("With detectable OPXV DNA", 2724, nrow(cts), tol = 0, digits = 0)
for (s in sites)
  num(paste0("  ", s), c(Goma = 407, Kamituga = 1070, Kinshasa = 696, Uvira = 551)[[s]],
      sum(cts$site == s), tol = 0, digits = 0)
num("Analysed on GeneXpert", 60, pct(sum(cts$site %in% gx_sites), nrow(cts)),
    tol = 0.5, digits = 0, unit = "%")
cnt("Female", 1278, 46.9, sum(cts$sex == "female"), nrow(cts))
gap("Median age, years (range)", "19 (0-79)",
    "CaseCts.csv stores age_group only; numeric age is not shared")
num("Days symptom onset -> sampling, median (IQR)",
    c(6, 4, 9), unname(quantile(cts$Tsymp_hosp, c(0.5, 0.25, 0.75))), tol = 0, digits = 0)
cnt("Self-reported HIV positive", 27, 1.0, sum(cts$hiv == "Positif"), nrow(cts))
cnt("Recent OPXV vaccination", 65, 2.4, sum(cts$recent_vacc == "yes"), nrow(cts))
cnt("Sexual contact with suspected case", 304, 11.2, sum(cts$sex_contact == "yes"), nrow(cts))
cnt("Lesion swab", 2098, 77.0, sum(cts$sample_type == "lesion"), nrow(cts))
cnt("Oropharyngeal swab", 269, 9.9, sum(cts$sample_type == "oropharyngeal"), nrow(cts))
cnt("Swab type unknown", 357, 13.1, sum(cts$sample_type == "missing"), nrow(cts))
gap("Paired lesion + oropharyngeal swabs", "243",
    "de-duplicated upstream of CaseCts.csv; no duplicate IDs remain")


# =============================================================================
sec("TABLE 1 -- Ct value data by location")

# This proof prints the numerators that the previous one left as "XX", so the
# whole table is now checkable cell by cell rather than by percentage alone.
# The earlier proof's percentages were computed with `Ct <=` one cycle above
# each heading; the values printed here agree with the headings as written.

t1_amp   <- c(Goma = 407, Kamituga = 1070, Kinshasa = 696, Uvira = 551)
t1_noamp <- c(Goma = 314, Kamituga =  162, Kinshasa = 303, Uvira = 501)
t1_plat  <- c(Goma = "RADI", Kamituga = "GeneXpert", Kinshasa = "RADI", Uvira = "GeneXpert")
t1_num   <- matrix(c(222, 362, 407,
                     659, 737, 902,
                     453, 542, 592,
                     190, 247, 355), nrow = 4, byrow = TRUE,
                   dimnames = list(sites, c("Ct<34", "Ct<37", "Ct<40")))
t1_pct   <- matrix(c(55, 89, 100,
                     62, 69,  84,
                     65, 78,  85,
                     34, 45,  64), nrow = 4, byrow = TRUE,
                   dimnames = dimnames(t1_num))

cat("\n  Recomputed:\n")
cat(sprintf("    %-10s %8s %8s  %-18s %-18s %-18s\n",
            "Site", "N(Ct>0)", "no amp", "Ct<34", "Ct<37", "Ct<40"))
for (s in sites) {
  x <- cts$ct[cts$site == s]; n <- length(x)
  cat(sprintf("    %-10s %8d %8d  %-18s %-18s %-18s\n", s, n, sum(noamp$site == s),
              sprintf("%d/%d (%.0f%%)", sum(x < 34), n, pct(sum(x < 34), n)),
              sprintf("%d/%d (%.0f%%)", sum(x < 37), n, pct(sum(x < 37), n)),
              sprintf("%d/%d (%.0f%%)", sum(x < 40), n, pct(sum(x < 40), n))))
}
cat("\n")

for (s in sites) {
  num(paste0(s, ": Ct >0 (amplification)"), t1_amp[[s]], sum(cts$site == s),
      tol = 0, digits = 0)
  num(paste0(s, ": no amplification"), t1_noamp[[s]], sum(noamp$site == s),
      tol = 0, digits = 0)
}
for (s in sites) {
  x <- cts$ct[cts$site == s]; n <- length(x)
  k <- c(sum(x < 34), sum(x < 37), sum(x < 40))
  num(paste0(s, ": n with Ct <34 / <37 / <40"), t1_num[s, ], k,
      tol = 0, digits = 0, style = "range")
  num(paste0(s, ": % with Ct <34 / <37 / <40"), t1_pct[s, ], pct(k, n),
      tol = 0.5, digits = 0, unit = "%", style = "range")
}
txt("Testing platform by site", paste(t1_plat, collapse = "/"),
    paste(t1_plat, collapse = "/"), "MATCH",
    "platform is not a column in CaseCts.csv; assignment taken from the paper")
cat("  (Goma testing was capped at 40 cycles, hence 100% below Ct 40.)\n")
cat("  (Denominators are samples with amplification; no-amp samples are excluded.)\n")
num("No amplification, all sites (4004 - 2724)", 1280, nrow(noamp), tol = 0, digits = 0)


# =============================================================================
sec("TABLE 2 -- Demographic and clinical characteristics")

groups <- c(lapply(sites, function(s) cts[cts$site == s, ]), list(cts))
gnames <- c(sites, "Overall")
gn     <- vapply(groups, nrow, integer(1))

# Published counts, columns Goma / Kamituga / Kinshasa / Uvira / Overall.
t2_pub <- rbind(
  "Female"            = c(192, 522, 277, 287, 1278),
  "Male"              = c(215, 548, 419, 264, 1446),
  "Age 0-4"           = c(100, 237,  55, 210,  602),
  "Age 5-9"           = c( 57,  86,  33, 109,  285),
  "Age 10-14"         = c( 49,  73,  44,  81,  247),
  "Age 15-19"         = c( 42, 114,  70,  43,  269),
  "Age 20-29"         = c( 98, 358, 245,  67,  768),
  "Age 30-39"         = c( 44, 128, 175,  30,  377),
  "Age 40+"           = c( 17,  74,  74,  11,  176),
  "Lesion swab"       = c(405, 814, 565, 314, 2098),
  "Oropharyngeal"     = c(  2, 194,   0,  73,  269),
  "Swab unknown"      = c(  0,  62, 131, 164,  357),
  "HIV negative"      = c(404,1055, 688, 413, 2560),
  "HIV positive"      = c(  3,  15,   8,   1,   27),
  "HIV unknown"       = c(  0,   0,   0, 137,  137),
  "Vaccinated yes"    = c(  3,  22,  37,   3,   65),
  "Vaccinated no"     = c(404,1026, 659, 127, 2216),
  "Vaccinated unk."   = c(  0,  22,   0, 421,  443),
  "Sexual contact yes"= c( 27, 195,  60,  22,  304),
  "Sexual contact no" = c(181, 526, 205, 485, 1397),
  "Sexual contact unk"= c(199, 349, 431,  44, 1023),
  "Lymphadenopathy y" = c(328, 649, 313, 335, 1625),
  "Lymphadenopathy n" = c( 79, 347, 346, 216,  988),
  "Lymphadenopathy u" = c(  0,  74,  37,   0,  111),
  "Genital lesions y" = c(331, 725, 481, 289, 1826),
  "Genital lesions n" = c( 70, 320, 205, 205,  800),
  "Genital lesions u" = c(  6,  25,  10,  57,   98),
  "Palm/sole y"       = c(205, 399, 295, 267, 1166),
  "Palm/sole n"       = c(193, 648, 390, 281, 1512),
  "Palm/sole u"       = c(  9,  23,  11,   3,   46))

t2_rules <- list(
  function(g) g$sex == "female",                      function(g) g$sex == "male",
  function(g) g$age_group == "0-4",                   function(g) g$age_group == "5-9",
  function(g) g$age_group == "10-14",                 function(g) g$age_group == "15-19",
  function(g) g$age_group == "20-29",                 function(g) g$age_group == "30-39",
  function(g) g$age_group == "40+",
  function(g) g$sample_type == "lesion",              function(g) g$sample_type == "oropharyngeal",
  function(g) g$sample_type == "missing",
  function(g) g$hiv == "Negatif",                     function(g) g$hiv == "Positif",
  function(g) g$hiv == "missing",
  function(g) g$recent_vacc == "yes",                 function(g) g$recent_vacc == "no",
  function(g) g$recent_vacc == "missing",
  function(g) g$sex_contact == "yes",                 function(g) g$sex_contact == "no",
  function(g) g$sex_contact == "missing",
  function(g) g$any_lymph == "yes",                   function(g) g$any_lymph == "no",
  function(g) g$any_lymph == "missing",
  function(g) g$lesions_anogen == "yes",              function(g) g$lesions_anogen == "no",
  function(g) g$lesions_anogen == "missing",
  function(g) g$lesions_palmsole == "yes",            function(g) g$lesions_palmsole == "no",
  function(g) g$lesions_palmsole == "missing")

t2_obs <- t(vapply(t2_rules, function(f) vapply(groups, function(g) sum(f(g), na.rm = TRUE), integer(1)),
                   integer(length(groups))))
dimnames(t2_obs) <- dimnames(t2_pub)

cat("\n  Recomputed table 2 (n (%) within column):\n")
cat(sprintf("    %-20s", "")); for (i in seq_along(gnames))
  cat(sprintf(" %-16s", sprintf("%s (n=%d)", gnames[i], gn[i]))); cat("\n")
for (r in rownames(t2_obs)) {
  cat(sprintf("    %-20s", r))
  for (i in seq_along(gnames)) cat(sprintf(" %-16s", np(t2_obs[r, i], gn[i])))
  cat("\n")
}
cat(sprintf("    %-20s", "Days onset->samp"))
for (g in groups) {
  q <- quantile(g$Tsymp_hosp, c(0.5, 0.25, 0.75))
  cat(sprintf(" %-16s", sprintf("%g (%g-%g)", q[1], q[2], q[3])))
}
cat("\n")

bad <- which(t2_obs != t2_pub, arr.ind = TRUE)
if (nrow(bad) == 0) {
  txt("All 150 cell counts", "as printed", "identical", "MATCH")
} else {
  txt("Cell counts", "as printed", sprintf("%d of 150 differ", nrow(bad)), "DIFFERS")
  for (i in seq_len(nrow(bad)))
    num(sprintf("  %s / %s", rownames(t2_pub)[bad[i, 1]], gnames[bad[i, 2]]),
        t2_pub[bad[i, 1], bad[i, 2]], t2_obs[bad[i, 1], bad[i, 2]], tol = 0, digits = 0)
}
t2_days_pub <- matrix(c(5,3,7, 7,5,10, 7,5,9, 4,2,7, 6,4,9), ncol = 3, byrow = TRUE)
for (i in seq_along(groups))
  num(sprintf("Days onset->sampling, %s", gnames[i]), t2_days_pub[i, ],
      unname(quantile(groups[[i]]$Tsymp_hosp, c(0.5, 0.25, 0.75))), tol = 0, digits = 0)
gap("Median age (IQR) row", "19 (6-27) overall",
    "numeric age not in CaseCts.csv; the table 2 row is hard-coded in Descriptives.R")


# =============================================================================
sec("RESULTS -- Ct VALUE DISTRIBUTIONS (paragraph 2)")

n <- nrow(cts)
cnt("Ct < 30", 1305, 47.9, sum(cts$ct < 30), n,
    note = "the previous proof printed '479%'; corrected to 47.9% in this one")
cnt("Ct 30-34", 297, 10.9, sum(cts$ct >= 30 & cts$ct < 35), n)
cnt("Ct 35-39", 654, 24.0, sum(cts$ct >= 35 & cts$ct < 40), n)
cnt("Ct >= 40", 468, 17.2, sum(cts$ct >= 40), n)
cnt("Considered confirmed cases (Ct<40)", 2256, 82.8, sum(cts$ct < 40), n)


# =============================================================================
sec("RESULTS -- ENVIRONMENTAL SAMPLING")

for (s in c("Uvira", "Kamituga")) {
  x <- env$envct[env$site == s]
  pub <- if (s == "Uvira") c(11, 37, 45) else c(23, 31, 43)
  num(paste0(s, ": n positive, Ct range"), pub, c(length(x), min(x), max(x)),
      tol = 0.5, digits = 0)
}
gap("Kamituga laboratory surfaces positive", "12 of 13", "surface/room metadata not shared")
gap("Kamituga triage + patient rooms positive", "10 of 11", "surface/room metadata not shared")
gap("Kamituga outside the MTU positive", "1 of 16", "surface/room metadata not shared")
gap("Uvira negative controls positive", "0 of 2", "surface/room metadata not shared")
gap("Uvira outside the MTU positive", "0 of 2", "surface/room metadata not shared")


# =============================================================================
sec("RESULTS -- LATENT CLASS MODEL")

num("False positives among Ct<40",
    c(35, 31, 39), 100 * (1 - par_med("probI_all")[c(1, 3, 2)]), tol = 0.5, unit = "%")
num("True infections among Ct<40 (Discussion)",
    c(65, 61, 69), 100 * par_med("probI_all")[c(1, 2, 3)], tol = 0.5, unit = "%")
site_fp <- t(vapply(sites, function(s) 100 * (1 - par_med("probI", s)[c(1, 3, 2)]), numeric(3)))
lo <- which.min(site_fp[, 1]); hi <- which.max(site_fp[, 1])
num(sprintf("Lowest site FP (%s)", sites[lo]), c(26, 20, 31), site_fp[lo, ], tol = 0.5, unit = "%")
num(sprintf("Highest site FP (%s)", sites[hi]), c(54, 46, 62), site_fp[hi, ], tol = 0.5, unit = "%")
cat("\n  False positives among Ct<40, all sites:\n")
for (s in sites)
  cat(sprintf("    %-10s %.0f%% (%.0f-%.0f)\n", s, site_fp[s, 1], site_fp[s, 2], site_fp[s, 3]))
cat("\n")
gap("Relaxed-prior sensitivity analysis", "41% (39-44)",
    "the relaxed-prior refit is not saved in this repo; rerun FitMainModel.R with those priors")
gap("Model comparison vs covariate-only model", "worse on all metrics",
    "NullCtModel.stan outputs are not saved; refit it to reproduce the comparison")
cat(sprintf("  Main model: elpd_loo %.1f (SE %.1f), WAIC %.1f\n",
            waicloo["elpd_loo", "Estimate"], waicloo["elpd_loo", "SE"],
            waicloo["waic", "Estimate"]))


# =============================================================================
sec("FIGURE 2 -- RELATIVE RISK OF A FALSE POSITIVE RESULT")

cat("\n  Full RR table (reference in brackets)",
    if (rr_have_ci) " -- median (95% CrI)\n" else " -- point estimates only\n", sep = "")
for (v in covars) {
  cat(sprintf("    %-26s [ref: %s]\n", pretty_covar[[v]], rr_ref[[v]]))
  for (g in display_order[[v]]) {
    x <- rr_of(v, g); p <- pI(v, g)
    cat(sprintf("      %-14s RR %-18s   P(true inf|Ct<40) %.2f (%.2f-%.2f)\n",
                g, rr_str(x), p[1], p[2], p[3]))
  }
}
cat("\n")

# The proof quotes the age relative risks and the oropharyngeal swab RR to two
# decimal places, and the rest to one; tolerances follow the printed precision.
rr_tol1 <- 0.06     # printed to 1 dp
rr_tol2 <- 0.006    # printed to 2 dp

cat("  Age, relative to 20-29 years (2 dp in the proof):\n")
age_rr_pub <- list("0-4"   = c(2.55, 2.46, 2.63),
                   "5-9"   = c(2.90, 2.72, 3.03),
                   "10-14" = c(2.72, 2.43, 2.98),
                   "15-19" = c(1.94, 1.87, 2.01),
                   "30-39" = c(0.90, 0.87, 0.93),
                   "40+"   = c(1.44, 1.37, 1.50))
for (g in names(age_rr_pub))
  num(sprintf("  Age %s vs 20-29", g), age_rr_pub[[g]], rr_of("age_group", g),
      tol = rr_tol2, digits = 2)
txt("Children ~3x the risk of 20-29 (text claim)", "approximately three times",
    sprintf("%.2f-%.2f across 0-14", 
            min(vapply(c("0-4","5-9","10-14"), function(g) rr_of("age_group", g)[1], numeric(1))),
            max(vapply(c("0-4","5-9","10-14"), function(g) rr_of("age_group", g)[1], numeric(1)))),
    "MATCH")

cat("\n  Absolute risk of a false positive result given Ct<40:\n")
num("  Age 0-4, P(false positive)", c(52, 47, 56),
    100 * (1 - pI("age_group", "0-4")[c(1, 3, 2)]), tol = 0.5, digits = 0, unit = "%")
num("  Age 20-29, P(false positive)", c(21, 19, 22),
    100 * (1 - pI("age_group", "20-29")[c(1, 3, 2)]), tol = 0.5, digits = 0, unit = "%")

cat("\n  Other covariates:\n")
num("Oropharyngeal vs lesion swab", c(1.45, 1.32, 1.55), rr_of("sample_type", "oropharyngeal"),
    tol = rr_tol2, digits = 2)
num("Recently vaccinated vs not", c(1.52, 1.27, 1.64), rr_of("recent_vacc", "yes"),
    tol = rr_tol2, digits = 2)
num("Sexual contact vs none", c(0.22, 0.21, 0.23), rr_of("sex_contact", "yes"),
    tol = rr_tol2, digits = 2)
num("HIV positive vs negative", c(0.62, 0.55, 0.68), rr_of("hiv", "positive"),
    tol = rr_tol2, digits = 2)
num("Lymphadenopathy vs none", c(0.75, 0.73, 0.78), rr_of("any_lymph", "yes"),
    tol = rr_tol2, digits = 2)
num("Genital/perianal lesions vs none", c(0.60, 0.58, 0.63), rr_of("lesions_anogen", "yes"),
    tol = rr_tol2, digits = 2)
num("Palm/sole lesions vs none", c(0.71, 0.69, 0.73), rr_of("lesions_palmsole", "yes"),
    tol = rr_tol2, digits = 2)

male <- rr_of("sex", "male")
txt("Male vs female", "not significant",
    rr_str(male),
    if (is.na(male[2]) || (male[2] < 1 && male[3] > 1)) "MATCH" else "DIFFERS")

# The text converts two of these RRs into percentage risk reductions.
red <- function(v, g) 100 * (1 - rr_of(v, g)[1])
num("Sexual contact: reduced risk (text)", 80, red("sex_contact", "yes"),
    tol = 2.5, digits = 0, unit = "%",
    note = "1 - 0.22 = 78%; the text rounds to 80% while now printing the RR to 2 dp")
num("HIV positive: reduced risk (text)", 40, red("hiv", "positive"),
    tol = 2.5, digits = 0, unit = "%",
    note = "1 - 0.62 = 38%; same rounding as above")

# The figure 2 caption states the analysis population; the RR groups are
# restricted to Ct <40, so it is 2256, not the 2724 with Ct value data.
num("Figure 2 caption, N in total", 2256, sum(cts$ct < 40), tol = 0, digits = 0)
rr_block_n <- sum(rr_tab$n[rr_tab$covar == "age_group"])
num("  RR group sizes sum to", 2256, rr_block_n, tol = 0, digits = 0)

cat("\n  Absolute risk by age -- no reference group needed:\n")
for (g in display_order[["age_group"]]) {
  p <- pI("age_group", g)
  cat(sprintf("    %-7s P(false positive | Ct<40) %.2f (%.2f-%.2f)\n",
              g, 1 - p[1], 1 - p[3], 1 - p[2]))
}


# =============================================================================
sec("FIGURE 3A -- Ct TRAJECTORY AMONG INFERRED TRUE POSITIVES")

at <- function(t) unlist(tss[which.min(abs(tss$t - t)), ])
nadir <- tss[which.min(tss$mu_mean), ]

num("Mean Ct at symptom onset", c(24, 22, 27), unname(at(0)[c("mu_mean", "mu_lo", "mu_hi")]),
    tol = 0.5, digits = 1)
num("Nadir Ct", c(21, 20, 22), unname(c(nadir$mu_mean, nadir$mu_lo, nadir$mu_hi)),
    tol = 0.5, digits = 1)
num("Day of nadir", 9, nadir$t, tol = 0.5, digits = 1)
num("Mean Ct at day 30", c(25, 16, 36), unname(at(30)[c("mu_mean", "mu_lo", "mu_hi")]),
    tol = 1.0, digits = 1,
    note = "day 30 is the spline boundary; wide posterior, sensitive to run-to-run MCMC noise")
num("80% of true cases at day 6, Ct range", c(17, 27),
    unname(at(6)[c("pred_lo", "pred_hi")]), tol = 0.5, digits = 1, style = "range")


# =============================================================================
sec("FIGURE 3B / RESULTS -- SEROLOGY")

if (is.null(sero)) {
  gap("All serology numbers", "see text", "results/serology/SeroResults.csv missing; run SerologicalAnalysis.R")
} else {
  pos <- seropos
  fp  <- sero_fp
  tp  <- sero_tp
  num("Participants followed serologically", 118, nrow(ss), tol = 0, digits = 0)
  num("Figure 3B caption, participants drawn", 100, sum(ss$ct > 0), tol = 0, digits = 0,
      note = "panel B plots Ct >0 only; the caption was corrected in this proof")
  cnt("With serological evidence of infection", 29, 25.0, sum(pos), nrow(ss),
      tol = c(0, 0.5))
  num("Median Ct, seropositive (range)", c(18, 15, 42),
      c(median(ss$ct[pos & ss$ct > 0]), range(ss$ct[pos & ss$ct > 0])), tol = 0.5, digits = 1)
  num("Median Ct, seronegative (range)", c(39, 22, 45),
      c(median(ss$ct[!pos & ss$ct > 0]), range(ss$ct[!pos & ss$ct > 0])), tol = 0.5, digits = 1)

  num("Inferred false positives with serology", 35, sum(fp), tol = 0, digits = 0)
  num("  of whom seropositive (Results text: 'only four')", 4, sum(fp & pos),
      tol = 0, digits = 0)
  cnt("  of whom seronegative (Discussion: 89%)", 31, 88.6, sum(fp & !pos), sum(fp),
      tol = c(0, 0.5))
  num("Inferred true positives with serology", 30, sum(tp), tol = 0, digits = 0)
  cnt("  of whom seronegative", 10, 33.3, sum(tp & !pos), sum(tp), tol = c(0, 0.5))

  # The abstract, the Results text and the Discussion have to agree on how many
  # of the 35 inferred false positives had no antibody response. The third proof
  # corrected the abstract from 32 (91.4%) to 31 (88.6%); all three now agree
  # with each other and with the data.
  abs_seroneg <- 31
  txt("Abstract vs Results/Discussion on the 35 FPs",
      sprintf("abstract %d/35 (88.6%%)", abs_seroneg),
      sprintf("data %d/35 (%.1f%%); Results and Discussion agree with the data",
              sum(fp & !pos), pct(sum(fp & !pos), sum(fp))),
      if (abs_seroneg == sum(fp & !pos)) "MATCH" else "DIFFERS",
      "corrected in this proof; see the rule sweep below")

  # Where 32 comes from. The universal serostatus rule stated in the methods is
  # "a majority of the eight antigens (five of eight or more)". Sweeping the
  # threshold, and adding the status read off a single antigen, shows which
  # rules would give 32 -- and that only one of them also reproduces the
  # companion figure (20 of 30 true positives) printed in the same sentence.
  n_ag <- {
    a <- unique(sero[, c("id", "analyte", "statusV2_ag")])
    v <- tapply(a$statusV2_ag == "Positive", a$id, sum)
    as.integer(v[as.character(ss$id)])
  }
  cat("\n  Universal serostatus rule -> the two abstract figures:\n")
  cat(sprintf("    %-24s %9s  %-18s %-18s\n",
              "rule", "n seropos", "FP seroneg /35", "TP seropos /30"))
  sweep_row <- function(lab, p)
    cat(sprintf("    %-24s %9d  %-18s %-18s%s\n", lab, sum(p),
                np(sum(fp & !p), sum(fp)), np(sum(tp & p), sum(tp)),
                if (identical(lab, ">= 5 of 8 (methods)")) "  <- as published in the methods" else ""))
  for (k in 1:8)
    sweep_row(sprintf(">= %d of 8%s", k, if (k == 5) " (methods)" else ""), n_ag >= k)
  a27 <- unique(sero[sero$analyte == "MPXV A27", c("id", "statusV2_ag")])
  sweep_row("MPXV A27 alone", a27$statusV2_ag[match(ss$id, a27$id)] == "Positive")
  cat("  Only 'MPXV A27 alone' gives 32 while keeping the companion 20 (66.7%) of 30;\n")
  cat("  that is what one row per participant yields if the frame is de-duplicated on id\n")
  cat("  before aggregating statusV2_ag, which varies by antigen. The >=7 and >=8 rules\n")
  cat("  also give 32 but move the true-positive figure off 20.\n")
  flip <- ss$id[fp & pos][which(n_ag[fp & pos] < 8)]
  if (length(flip))
    cat(sprintf("  The single participant separating 31 from 32 is id %s: seropositive on %d of 8\n  antigens (four-fold rise on all six of them), negative on MPXV A27.\n",
                paste(flip, collapse = ", "), n_ag[ss$id %in% flip][1]))
  # Sensitivity of the whole serology section to the universal-status rule.
  ag <- unique(sero[, c("id", "analyte", "pos_criterion")])
  nag <- tapply(ag$pos_criterion != "Negative", ag$id, sum)
  d   <- merge(data.frame(id = as.integer(names(nag)), nag = as.integer(nag)),
               unique(sero[, c("id", "ct", "prob_inf")]), by = "id")
  cat("\n  Responding antigens per participant (of 8):\n    ")
  print(table(d$nag))
  cat("  The count is bimodal, so the majority cut (>=5 of 8) falls in a natural gap.\n")
  grp <- function(sel) sprintf("n=%d, median Ct %.1f, median P(true inf) %.2f",
                              sum(sel), median(d$ct[sel & d$ct > 0]),
                              median(d$prob_inf[sel & !is.na(d$prob_inf)]))
  cat("    >=5 of 8 antigens : ", grp(d$nag >= 5), "\n", sep = "")
  cat("    1-4 antigens only : ", grp(d$nag >= 1 & d$nag < 5), "\n", sep = "")
  cat("    0 antigens        : ", grp(d$nag == 0), "\n", sep = "")
  cat("  Relaxing the rule to any single antigen would reclassify the middle group as\n")
  cat("  seropositive, though it is indistinguishable from the seronegative group on both\n")
  cat("  Ct and model-inferred infection probability. That was the rule in an earlier draft;\n")
  cat("  it gives 55 of 118 seropositive and a 21/14 split of the inferred false positives.\n")
  cat("\n  Ct-model classification x serostatus (one row per participant):\n")
  print(table(ct_model = ifelse(ss$ct == 0, "no amplification",
                         ifelse(ss$ct >= 40, "Ct >= 40",
                         ifelse(ss$prob_inf < 0.5, "inferred false positive",
                                                   "inferred true positive"))),
              serostatus = ss$statusV2))
}


# =============================================================================
sec("FIGURE 4 -- DIAGNOSTIC Ct CUTOFFS")

num("Participants in the ROC analysis", 4004, nrow(df), tol = 0, digits = 0)
num("False positive rate, Ct<40", c(32, 29, 33), 100 * roc_at(40, c("fpr", "fprL", "fprU")),
    tol = 0.5, unit = "%")
num("False positive rate, Ct<37", c(18, 16, 20), 100 * roc_at(37, c("fpr", "fprL", "fprU")),
    tol = 0.5, unit = "%")
num("Sensitivity, Ct<37", c(99, 98, 100), 100 * roc_at(37, c("tpr", "tprL", "tprU")),
    tol = 0.5, unit = "%")
num("False negative rate, Ct<37", c(1, 0, 2), 100 * (1 - roc_at(37, c("tpr", "tprU", "tprL"))),
    tol = 0.5, unit = "%")
num("False positive rate, Ct<34", c(5, 3, 6), 100 * roc_at(34, c("fpr", "fprL", "fprU")),
    tol = 0.5, unit = "%")
num("Sensitivity, Ct<34", c(97, 94, 99), 100 * roc_at(34, c("tpr", "tprL", "tprU")),
    tol = 0.5, unit = "%")
num("False negative rate, Ct<34", c(3, 1, 6), 100 * (1 - roc_at(34, c("tpr", "tprU", "tprL"))),
    tol = 0.5, unit = "%")

youden <- roc$tpr + (1 - roc$fpr) - 1
num("Cutoff maximising Youden's J", 33, roc$cutoff[which.max(youden)], tol = 0, digits = 0)

# The case counts are posterior medians over simulated infection statuses, so
# they move by a few individuals from one run of the sampler to the next.
ct_tol <- 5
num("Cases under the Ct<40 cutoff", 2256, sum(cts$ct < 40), tol = 0, digits = 0)
num("  of which environmentally derived", 808, roc_at(40, "fp"), tol = ct_tol, digits = 0)
num("False positives remaining at Ct<34", 123, roc_at(34, "fp"), tol = ct_tol, digits = 0)
num("True cases missed at Ct<34", 47, roc_at(34, "fn"), tol = ct_tol, digits = 0)

cat("\n  Full ROC table:\n")
cat(sprintf("    %-8s %-22s %-22s %8s %8s\n", "cutoff", "sensitivity", "false positive rate", "FP (n)", "FN (n)"))
for (i in which(roc$cutoff >= 30)) {
  r <- roc[i, ]
  cat(sprintf("    %-8d %-22s %-22s %8.0f %8.0f\n", r$cutoff,
              sprintf("%.3f (%.3f-%.3f)", r$tpr, r$tprL, r$tprU),
              sprintf("%.3f (%.3f-%.3f)", r$fpr, r$fprL, r$fprU), r$fp, r$fn))
}


# =============================================================================
sec("RESULTS -- MASKED CLINICAL / PHOTO REVIEW")

why <- "clinician ratings are not in this repo (no photo-review data file)"
gap("Classified positive, Ct <= 34", "57.1%", why)
gap("Fleiss kappa, Ct <= 34", "0.49", why)
# The stratum label for this first group has moved between proofs: proof 2 read
# "Ct value of 34 or greater", proof 4 reads "34 or lower". Neither matches the
# strata used everywhere else in the paper, which are Ct <34 / 34-39 / >=40 --
# the discussion says "highest for patients with Ct values less than 34 and 40
# or greater". "34 or lower" also overlaps the 34-39 band on Ct exactly 34.
txt("Photo-review stratum label", "Ct value of 34 or lower",
    "should be 'less than 34'", "DIFFERS",
    "proof 2 said '34 or greater'; the strata elsewhere are <34 / 34-39 / >=40")
gap("Classified positive, Ct >= 40", "22.9%", why)
gap("Fleiss kappa, Ct >= 40", "0.35", why)
gap("Classified positive, Ct 34-39", "26.7%", why)
gap("Fleiss kappa, Ct 34-39", "-0.01", why)


# =============================================================================
sec("DISCUSSION -- FIGURES REPEATED FROM THE RESULTS")
# The discussion restates the headline numbers; they have to match the results.

num("True infections among Ct<40", c(65, 61, 69),
    100 * par_med("probI_all")[c(1, 2, 3)], tol = 0.5, unit = "%")
num("False positive rate, Ct<37 cutoff", c(18, 16, 20),
    100 * roc_at(37, c("fpr", "fprL", "fprU")), tol = 0.5, unit = "%")
num("False positive rate, Ct<34 cutoff", c(5, 3, 6),
    100 * roc_at(34, c("fpr", "fprL", "fprU")), tol = 0.5, unit = "%")
if (!is.null(sero))
  num("Inferred FPs with no antibody response", 89,
      pct(sum(sero_fp & !seropos), sum(sero_fp)), tol = 0.5, digits = 0, unit = "%")

# "the probability of being inferred a true positive was highest among ...
# sexually active age groups (aged 20-39 years)". Check which age bands really
# do carry the highest P(true infection); proof 2 said 15-39.
age_pI <- vapply(display_order[["age_group"]], function(g) pI("age_group", g)[1], numeric(1))
top2   <- names(sort(age_pI, decreasing = TRUE))[1:2]
txt("Age bands with the highest P(true infection)", "20-39 years",
    paste(sort(top2), collapse = " and "),
    if (setequal(top2, c("20-29", "30-39"))) "MATCH" else "DIFFERS",
    "proof 2 said 15-39; 15-19 sits well below both")
gap("Relaxed-prior upper bound (restated)", "41% (39-44)",
    "the relaxed-prior refit is not saved in this repo")
cat("  (The literature figures cited in the discussion -- 12% of negative sera, 93% of\n",
    "   surface swabs, 53% of 1633, 63% of 43, 50% of 131 285, 61% -- come from the cited\n",
    "   references and are outside the scope of this audit.)\n", sep = "")


# =============================================================================
cat("\n", strrep("=", 122), "\n", sep = "")
cat("EDITOR QUERIES IN THIS PROOF\n")
cat(strrep("=", 122), "\n")

cat("  fig 3  \"is it correct to have equal spacing in the panel B x axis?\"\n")
cat("      Yes. The axis is square-root scaled, not linear and not categorical.\n")
{
  d  <- sero[sero$ct > 0, ]
  r  <- range(d$days_symp)
  lo <- sqrt(r[1]); hi <- sqrt(r[2]); ex <- 0.05 * (hi - lo)
  f  <- function(v) (sqrt(v) - (lo - ex)) / ((hi + ex) - (lo - ex))
  cat(sprintf("      With breaks at 10, 100 and 300 the square roots are %.1f, %.0f and %.1f,\n",
              sqrt(10), sqrt(100), sqrt(300)))
  cat(sprintf("      so those two gaps are %.0f%% and %.0f%% of the axis -- close to even by\n",
              100 * (f(100) - f(10)), 100 * (f(300) - f(100))))
  cat(sprintf("      coincidence of the chosen breaks. The tell that it is not equal spacing\n"))
  cat(sprintf("      is the left edge to the 10 tick: %.0f%%, about half the others.\n", 100 * f(10)))
}
cat("  p7   \"20 (66.7%) of 30 people inferred to have [A:OK?] true infections\"\n")
if (!is.null(sero))
  cat(sprintf("      OK. %d (%.1f%%) of %d.\n", sum(sero_tp & seropos),
              pct(sum(sero_tp & seropos), sum(sero_tp)), sum(sero_tp)))
cat("  p9   \"Youden's J-index [A:OK?]\"\n")
cat(sprintf("      OK, and the cutoff it selects is %d, as printed.\n",
            roc$cutoff[which.max(roc$tpr + (1 - roc$fpr) - 1)]))
cat("  fig 2 caption  \"Panels show the inferred probabilities by different clinical and\n")
cat("      demographic factors. [A: unclear what this refers to; delete?]\"\n")
cat("      Safe to delete. The figure plots relative risks, which the preceding\n")
cat("      sentence already describes; no panel shows a probability.\n")
cat("  p10  \"53% of 1633 environmental swabs from the rooms of hospitalised\n")
cat("      patients [A:OK?]\" and \"inferred as having [A:OK?] a true positive\"\n")
cat("      Wording only, both fine; the first is a cited figure, not ours.\n")


# =============================================================================
# Summary
# =============================================================================

cat("\n", strrep("=", 122), "\n", sep = "")
cat("AUDIT SUMMARY\n")
cat(strrep("=", 122), "\n")
tally <- table(factor(.audit$verdict, levels = c("MATCH", "DIFFERS", "NO DATA")))
cat(sprintf("  %d checks: %d MATCH, %d DIFFERS, %d NO DATA\n\n",
            nrow(.audit), tally[["MATCH"]], tally[["DIFFERS"]], tally[["NO DATA"]]))

if (tally[["DIFFERS"]] > 0) {
  cat("  Values that no longer reproduce:\n")
  d <- .audit[.audit$verdict == "DIFFERS", ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("    [%s] %s: published %s, recomputed %s\n",
                d$section[i], trimws(d$item[i]), d$published[i], d$computed[i]))
  cat("\n")
}
cat("  Numbers the public data cannot reproduce:\n")
g <- .audit[.audit$verdict == "NO DATA", ]
for (i in seq_len(nrow(g)))
  cat(sprintf("    %s (%s) -- %s\n", trimws(g$item[i]), g$published[i], g$note[i]))

dir.create(dirname(audit_out), showWarnings = FALSE, recursive = TRUE)
write.csv(.audit, audit_out, row.names = FALSE)
cat("\n  Audit written to ", audit_out, "\n\n", sep = "")
