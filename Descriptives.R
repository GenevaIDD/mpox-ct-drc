# =============================================================================
# Descriptive statistics for results paragraph
# =============================================================================
# Generates all in-text numbers for the opening Results paragraph.
#
# GAPS (data not available in CaseCts.csv):
#   - Numeric age: only age_group (categorical) is stored. The median age (IQR)
#     row in Table 2 is hard-coded from the upstream patient-level dataset,
#     which is not shared publicly; it cannot be recomputed from this repo.
#   - N=243 paired lesion+OP swabs: de-duplication was applied before saving
#     CaseCts.csv (no duplicate IDs), so this count is not recoverable here.

repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")

source(file.path(repo_dir, "LoadCaseCts.R"))

df  <- read_case_cts(data_dir)
cts <- df[df$ct > 0, ]          # Ct-positive (Ct > 0)
noamp <- df[df$ct == 0, ]       # no amplification

sites <- c("Goma", "Kamituga", "Kinshasa", "Uvira")

cat("========================================================\n")
cat("DESCRIPTIVE STATISTICS — Results paragraph\n")
cat("========================================================\n\n")

cat("--- Sample counts ---\n")
cat("Total qPCR results (Ct>0 and no-amp):", nrow(df), "\n")
cat("With detectable OPXV DNA (Ct>0):", nrow(cts), "\n")
cat("No amplification (Ct=0):", nrow(noamp), "\n")
cat("\nSite breakdown (Ct>0):\n")
print(table(cts$site)[sites])

cat("\n--- Testing platform ---\n")
n_genexpert <- sum(cts$site %in% c("Kamituga", "Uvira"))
cat("GeneXpert (Kamituga + Uvira):", n_genexpert,
    sprintf("(%d%%)\n", round(n_genexpert / nrow(cts) * 100)))
n_radi <- sum(cts$site %in% c("Goma", "Kinshasa"))
cat("RADI (Goma + Kinshasa):", n_radi,
    sprintf("(%d%%)\n", round(n_radi / nrow(cts) * 100)))

cat("\n--- Demographics ---\n")
cat("Female:", sum(cts$sex == "female"),
    sprintf("(%d%%)\n", round(mean(cts$sex == "female") * 100)))
cat("NOTE: numeric age not in CaseCts.csv — see hard-coded median age (IQR) in Table 2\n")
cat("Age group distribution:\n")
print(table(cts$age_group))

cat("\n--- Time from symptom onset to sample collection (Tsymp_hosp) ---\n")
cat("Median:", median(cts$Tsymp_hosp), "days\n")
cat("IQR:", paste(quantile(cts$Tsymp_hosp, c(0.25, 0.75)), collapse = "-"), "\n")
cat("Range:", paste(range(cts$Tsymp_hosp), collapse = "-"), "\n")

cat("\n--- Clinical covariates ---\n")
cat("HIV positive (Positif):", sum(cts$hiv == "Positif"),
    sprintf("(%d%%)\n", round(mean(cts$hiv == "Positif") * 100)))
cat("Recent vaccination (yes):", sum(cts$recent_vacc == "yes"),
    sprintf("(%d%%)\n", round(mean(cts$recent_vacc == "yes") * 100)))
cat("Sexual contact with suspected case (yes):", sum(cts$sex_contact == "yes"),
    sprintf("(%d%%)\n", round(mean(cts$sex_contact == "yes") * 100)))

cat("\n--- Swab type ---\n")
cat("Lesion swab:", sum(cts$sample_type == "lesion"),
    sprintf("(%d%%)\n", round(mean(cts$sample_type == "lesion") * 100)))
cat("Oropharyngeal swab:", sum(cts$sample_type == "oropharyngeal"),
    sprintf("(%d%%)\n", round(mean(cts$sample_type == "oropharyngeal") * 100)))
cat("Missing swab type:", sum(cts$sample_type == "missing"),
    sprintf("(%d%%)\n", round(mean(cts$sample_type == "missing") * 100)))
cat("NOTE: N=243 with both lesion+OP swabs not directly recoverable from CaseCts.csv\n")
cat("      (de-duplication was applied upstream; no duplicate IDs exist)\n")

cat("\n--- Ct value distributions (among Ct>0, N =", nrow(cts), ") ---\n")
cat("Ct < 30:", sum(cts$ct < 30),
    sprintf("(%d%%)\n", round(mean(cts$ct < 30) * 100)))
cat("Ct 30-34 (>=30 and <35):", sum(cts$ct >= 30 & cts$ct < 35),
    sprintf("(%d%%)\n", round(mean(cts$ct >= 30 & cts$ct < 35) * 100)))
cat("Ct 35-39 (>=35 and <40):", sum(cts$ct >= 35 & cts$ct < 40),
    sprintf("(%d%%)\n", round(mean(cts$ct >= 35 & cts$ct < 40) * 100)))
cat("Ct >= 40:", sum(cts$ct >= 40),
    sprintf("(%d%%)\n", round(mean(cts$ct >= 40) * 100)))
cat("Ct <  40 (considered confirmed positive):", sum(cts$ct < 40),
    sprintf("(%d%%)\n", round(mean(cts$ct < 40) * 100)))

cat("\n--- Table 1: Ct percentages by site ---\n")
cat(sprintf("%-12s  %6s  %10s  %-16s  %-16s  %-16s\n", "Site", "N (Ct>0)", "No amp", "Ct<34", "Ct<37", "Ct<40"))
for (s in sites) {
  x    <- cts$ct[cts$site == s]
  n    <- length(x)
  namp <- sum(noamp$site == s)
  fmt  <- function(k) sprintf("%d/%d (%d%%)", k, n, round(k / n * 100))
  cat(sprintf("%-12s  %6d  %10d  %-16s  %-16s  %-16s\n",
              s, n, namp,
              fmt(sum(x < 34)),
              fmt(sum(x < 37)),
              fmt(sum(x < 40))))
}
cat("NOTE: Goma testing was limited to 40 cycles (no Ct>=40 in data)\n")

cat("\n========================================================\n")
cat("TABLE 2: Demographic and clinical characteristics\n")
cat("         Participants with detectable MPXV DNA (Ct>0)\n")
cat("========================================================\n\n")

# ===== Helpers =====
pct1 <- function(k, n) sprintf("%d (%s%%)", k, formatC(k / n * 100, digits = 1, format = "f"))
med_iqr <- function(x) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf("%g (%g–%g)", q[1], q[2], q[3])
}

# ===== Column groups: each site, then overall =====
cols   <- c(sites, "Overall")
groups <- lapply(cols, function(s) if (s == "Overall") cts else cts[cts$site == s, ])
ns     <- sapply(groups, nrow)
header_line <- function() {
  cat(sprintf("%-35s", ""))
  for (i in seq_along(cols)) cat(sprintf("  %-18s", sprintf("%s (n=%d)", cols[i], ns[i])))
  cat("\n")
}

row_n <- function(label, k_vec) {
  cat(sprintf("  %-33s", label))
  for (i in seq_along(cols)) cat(sprintf("  %-18s", pct1(k_vec[i], ns[i])))
  cat("\n")
}
row_med <- function(label, vals_list) {
  cat(sprintf("  %-33s", label))
  for (v in vals_list) cat(sprintf("  %-18s", med_iqr(v)))
  cat("\n")
}
# Prints pre-formatted strings (used for the hard-coded age row below).
row_str <- function(label, strs) {
  cat(sprintf("  %-33s", label))
  for (s in strs) cat(sprintf("  %-18s", s))
  cat("\n")
}
section <- function(label) cat(sprintf("%-35s\n", label))

# ===== Print table =====
header_line()

section("Sex")
row_n("Female", sapply(groups, function(g) sum(g$sex == "female")))
row_n("Male",   sapply(groups, function(g) sum(g$sex == "male")))

# Numeric age is deliberately excluded from the public CaseCts.csv (only
# age_group is shared), so these values are HARD-CODED from the upstream
# patient-level dataset rather than recomputed. They cover participants with
# Ct>0 and were produced with the same quantile() defaults as med_iqr() above.
# Kinshasa's 25.5 is the median of an even-sized group (n=696): the two
# middle ages are 25 and 26.
age_med_iqr <- c(Goma     = "14 (5–25)",
                 Kamituga = "20 (6–27)",
                 Kinshasa = "25.5 (18–33)",
                 Uvira    = "7 (2–16)",
                 Overall  = "19 (6–27)")

section("Age, years")
row_str("Median (IQR)", age_med_iqr[cols])

section("Age group, years")
age_levels <- AGE_GROUP_LEVELS
age_labels <- c("0–4", "5–9", "10–14", "15–19", "20–29", "30–39", "≥40")
for (j in seq_along(age_levels)) {
  row_n(age_labels[j], sapply(groups, function(g) sum(g$age_group == age_levels[j], na.rm = TRUE)))
}

section("Median days from symptom onset to sampling")
row_med("Median (IQR)", lapply(groups, function(g) g$Tsymp_hosp))

section("Sample type")
row_n("Lesion",         sapply(groups, function(g) sum(g$sample_type == "lesion")))
row_n("Oropharyngeal",  sapply(groups, function(g) sum(g$sample_type == "oropharyngeal")))
row_n("Unknown",        sapply(groups, function(g) sum(g$sample_type == "missing")))

section("HIV status")
row_n("Negative", sapply(groups, function(g) sum(g$hiv == "Negatif")))
row_n("Positive", sapply(groups, function(g) sum(g$hiv == "Positif")))
row_n("Unknown",  sapply(groups, function(g) sum(g$hiv == "missing")))

section("Recent orthopoxvirus vaccination")
row_n("Yes",     sapply(groups, function(g) sum(g$recent_vacc == "yes")))
row_n("No",      sapply(groups, function(g) sum(g$recent_vacc == "no")))
row_n("Unknown", sapply(groups, function(g) sum(g$recent_vacc == "missing")))

section("Sexual contact")
row_n("Yes",     sapply(groups, function(g) sum(g$sex_contact == "yes")))
row_n("No",      sapply(groups, function(g) sum(g$sex_contact == "no")))
row_n("Unknown", sapply(groups, function(g) sum(g$sex_contact == "missing")))

section("Lymphadenopathy")
row_n("Yes",     sapply(groups, function(g) sum(g$any_lymph == "yes")))
row_n("No",      sapply(groups, function(g) sum(g$any_lymph == "no")))
row_n("Unknown", sapply(groups, function(g) sum(g$any_lymph == "missing")))

section("Genital or perianal lesions")
row_n("Yes",     sapply(groups, function(g) sum(g$lesions_anogen == "yes")))
row_n("No",      sapply(groups, function(g) sum(g$lesions_anogen == "no")))
row_n("Unknown", sapply(groups, function(g) sum(g$lesions_anogen == "missing")))

section("Palm or sole lesions")
row_n("Yes",     sapply(groups, function(g) sum(g$lesions_palmsole == "yes")))
row_n("No",      sapply(groups, function(g) sum(g$lesions_palmsole == "no")))
row_n("Unknown", sapply(groups, function(g) sum(g$lesions_palmsole == "missing")))
