# =============================================================================
# Ct mpox manuscript figures
# =============================================================================
# Reads the model outputs written by FitMainModel.R (results/main) and the
# serology results written by SerologicalAnalysis.R, and writes every
# manuscript figure to results/figures.

# ===== Libraries =====
library(ggplot2)
library(cowplot)
library(sf)
library(ggrepel)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# ===== Paths =====
repo_dir    <- getwd()
data_dir    <- file.path(repo_dir, "data")
res_dir     <- file.path(repo_dir, "results/main")
serores_dir <- file.path(repo_dir, "results/serology")
fig_dir     <- file.path(repo_dir, "results/figures")
tab_dir     <- file.path(repo_dir, "results/tables")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

source(file.path(repo_dir, "LoadCaseCts.R"))

# ===== Data =====
df  <- read_case_cts(data_dir)
env <- read.csv(file.path(data_dir, "EnvironmentalSamplingCts.csv"))

noamp <- df[df$ct == 0, ]   # qPCR results with no amplification
cts   <- df[df$ct > 0, ]    # qPCR results with Ct values

sites <- levels(factor(cts$site))

# ===== Model results =====
pars  <- read.csv(file.path(res_dir, "Params.csv"))
pic_rr   <- read.csv(file.path(res_dir, "pI_covars.csv"))
tsi   <- read.csv(file.path(res_dir, "TimeSinceSymps.csv"))
alpha <- read.csv(file.path(res_dir, "covar_coef.csv"))
pii   <- read.csv(file.path(res_dir, "pii.csv"))
ROC   <- readRDS(file.path(res_dir, "ROC.RDS"))
sero <- read.csv(file.path(serores_dir, "SeroResults.csv"))

# ===== Shared theme =====
theme_mpox <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95"),
    strip.text        = element_text(size = 11)
  )

covar_levels <- c("Age", "Sex", "Sample type", "Sexual contact", "HIV",
                  "Vaccinated", "Lymphadenopathy",
                  "Genital/perianal lesions", "Palm/sole lesions")

group_levels <- c(
  "0-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40+",
  "female", "male",
  "lesion", "oropharyngeal",
  "no", "yes",
  "negative", "positive"
)


# =============================================================================
# Figure 1: Study overview — Ct distributions + map
# =============================================================================

# 1A: Ct histograms by site
f1a <- ggplot(cts, aes(ct)) +
  geom_histogram(bins = 35, fill = "yellowgreen", col="grey20", linewidth=0.3) +
  geom_vline(xintercept = 40, col = "steelblue", linetype = "dashed",
             linewidth = 0.6) +
  geom_rug(data = env, aes(envct), sides = "b", col = "tomato",
           alpha = 0.9, length = unit(0.06, "npc")) +
  facet_wrap(~site) +
  theme_mpox + coord_cartesian(clip = "off") +
  labs(x = "Ct value", y = "Count")

f1a

# 1B: DRC map
sites_df <- tibble::tibble(
  site = c("Kinshasa", "Goma", "Uvira", "Kamituga"),
  lon  = c(15.2663,   29.2350, 29.1410, 27.9960),
  lat  = c(-4.4419,   -1.6790, -3.3950, -3.0730)
)
sites_sf <- st_as_sf(sites_df, coords = c("lon", "lat"), crs = 4326)
drc      <- ne_countries(scale = "medium",
                         country = "Democratic Republic of the Congo",
                         returnclass = "sf")
bbox  <- st_bbox(drc)
pad_x <- (bbox$xmax - bbox$xmin) * 0.08
pad_y <- (bbox$ymax - bbox$ymin) * 0.08

f1b <- ggplot() +
  geom_sf(data = drc, fill = "grey95", color = "grey40", linewidth = 0.4) +
  geom_sf(data = sites_sf, shape = 21, size = 3.5, stroke = 0.7,
          fill = "purple") +
  geom_text_repel(
    data = sites_df, aes(x = lon, y = lat, label = site),
    size = 3.5, box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0, segment.color = "grey40"
  ) +
  coord_sf(xlim = c(bbox$xmin - pad_x, bbox$xmax + pad_x),
           ylim = c(bbox$ymin - pad_y, bbox$ymax + pad_y),
           expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    plot.margin  = margin(8, 8, 8, 8),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4)
  )

fig1 <- plot_grid(f1b, f1a, rel_widths = c(0.55, 1), align = "b",
                  labels = c("A", "B"))
fig1

# ===== Output =====
ggsave(file.path(fig_dir, "Fig1.pdf"), fig1,
       width = 18, height = 8, units = "cm", dpi = 300)


# =============================================================================
# Figure 2: Relative risk of false positive (computed from chains)
# =============================================================================
library(posterior)
library(cmdstanr)

# Reference groups (RR = 1). Change any entry to use a different reference.
ref_groups <- c(
  "Age"                      = "20-29",
  "Sex"                      = "female",
  "Sample type"              = "lesion",
  "Sexual contact"           = "no",
  "HIV"                      = "negative",
  "Vaccinated"               = "no",
  "Lymphadenopathy"          = "no",
  "Genital/perianal lesions" = "no",
  "Palm/sole lesions"        = "no"
)

# Load saved fit and data objects
fit_obj  <- readRDS(file.path(res_dir, "Fit.RDS"))
data_rds <- readRDS(file.path(res_dir, "Data.RDS"))

# Covariate column order (must match order used during fitting)
covars_model <- c("sex","sample_type","age_group","sex_contact","recent_vacc","hiv",
                  "lesions_anogen","any_lymph","lesions_palmsole")
covar_name_map <- c(
  sex="Sex", sample_type="Sample type", age_group="Age",
  sex_contact="Sexual contact", recent_vacc="Vaccinated", hiv="HIV",
  lesions_anogen="Genital/perianal lesions", any_lymph="Lymphadenopathy",
  lesions_palmsole="Palm/sole lesions"
)

# Per-individual P(true infection) chains from the Stan fit
wb40     <- which(exp(data_rds$y) < if (!is.null(data_rds$ct_threshold)) data_rds$ct_threshold else 40)
pc_draws <- draws_of(as_draws_rvars(fit_obj$draws("pC"))$pC)
pc_draws <- exp(pc_draws)
pIi      <- pc_draws[,,1] / apply(pc_draws, c(1, 2), sum)
pIi_sub  <- pIi[, wb40]   # S x N_sub: Ct<40 individuals only

# For each covariate, match chain groups to labels sequentially:
# both pic_rr rows and the g=1..NG[c] loop skip groups with <2 obs in the
# same order, so their non-skipped entries correspond 1:1.
rr_rows <- list()
for (cv in covars_model) {
  c_idx     <- which(covars_model == cv)
  cov_col   <- data_rds$cov[wb40, c_idx]
  covar_lab <- covar_name_map[cv]
  ref_lab   <- ref_groups[covar_lab]

  pic_cov   <- pic_rr[as.character(pic_rr$covar) == covar_lab, ]
  if (nrow(pic_cov) == 0) next

  fp_chains <- list()
  n_per_grp <- list()
  pic_row_i <- 0L
  for (g in seq_len(data_rds$NG[c_idx])) {
    idx <- which(cov_col == g)
    if (length(idx) < 2) next
    pic_row_i <- pic_row_i + 1L
    label <- as.character(pic_cov$var[pic_row_i])
    if (label == "missing") next
    fp_chains[[label]] <- 1 - rowMeans(pIi_sub[, idx, drop = FALSE])
    n_per_grp[[label]] <- length(idx)
  }

  if (!ref_lab %in% names(fp_chains)) next
  fp_ref <- fp_chains[[ref_lab]]

  for (grp in names(fp_chains)) {
    if (grp == ref_lab) {
      rr_med <- 1; rr_lo <- NA_real_; rr_hi <- NA_real_
    } else {
      q      <- quantile(fp_chains[[grp]] / fp_ref, c(0.5, 0.025, 0.975), na.rm = TRUE)
      rr_med <- q[1]; rr_lo <- q[2]; rr_hi <- q[3]
    }
    # Group false-positive probability itself (the numerator of the RR), kept
    # for the supplementary table so the ratios can be read against a scale.
    fp_q <- quantile(fp_chains[[grp]], c(0.5, 0.025, 0.975), na.rm = TRUE)
    rr_rows[[paste(cv, grp, sep = "|||")]] <- data.frame(
      covar = covar_lab, var = grp, ref = ref_lab,
      rr = rr_med, rr_lo = rr_lo, rr_hi = rr_hi,
      fp = fp_q[1], fp_lo = fp_q[2], fp_hi = fp_q[3],
      n = n_per_grp[[grp]]
    )
  }
}

d_rr <- do.call(rbind, rr_rows)
d_rr$covar <- factor(d_rr$covar, levels = covar_levels)
d_rr$var   <- factor(d_rr$var,   levels = intersect(group_levels, unique(d_rr$var)))
d_rr       <- d_rr[!is.na(d_rr$covar), ]

# Every group, reference rows included: the supplementary table below reports
# these, while the figure drops the reference row of the binary covariates.
d_rr_full  <- d_rr[order(d_rr$covar, d_rr$var), ]
rownames(d_rr_full) <- NULL

# Drop reference row for binary covariates
n_grps      <- tapply(as.character(d_rr$var), as.character(d_rr$covar),
                      function(x) length(unique(x)))
binary_covs <- names(n_grps)[n_grps == 2]
ref_keys    <- paste(names(ref_groups), ref_groups, sep = "|||")
d_rr$key    <- paste(as.character(d_rr$covar), as.character(d_rr$var), sep = "|||")
d_rr        <- d_rr[!(d_rr$key %in% ref_keys & as.character(d_rr$covar) %in% binary_covs), ]
d_rr$key    <- NULL

# shared layers for fig 2 variants
fig2_base <- list(
  geom_vline(xintercept = 1, linetype = "dashed", col = "grey30", linewidth = 0.5),
  geom_linerange(aes(xmin = rr_lo, xmax = rr_hi), col = "blue", linewidth = 0.5),
  facet_grid(covar ~ ., scales = "free_y", space = "free", switch = "y"),
  scale_x_continuous(trans  = "log2",
                     breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3),
                     labels = c("0.25","0.5","0.75","1","1.5","2","3")),
  theme_mpox,
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10),
        strip.placement    = "outside",
        axis.title.y       = element_blank(),
        panel.spacing.y    = unit(2, "pt")),
  xlab("RR false positive diagnosis"),
  scale_y_discrete(position = "right")
)

# ===== Fig 2: original (no sample size) =====
fig2_rr <- ggplot(d_rr, aes(y = var, x = rr)) +
  fig2_base +
  geom_point(col = "blue", size = 2)
fig2_rr

# ===== Fig 2 variant A: dot size proportional to n =====
fig2_v_size <- ggplot(d_rr, aes(y = var, x = rr)) +
  fig2_base +
  geom_point(aes(size = n), col = "blue") +
  scale_size_continuous(name = "n", range = c(1, 5))
fig2_v_size


# output
ggsave(file.path(fig_dir, "Fig2.pdf"), fig2_v_size,
       width = 6.5, height = 7, units = "in", dpi = 300)


# =============================================================================
# Supplementary table: the numbers behind Figure 2
# =============================================================================
# One row per covariate group among the Ct<40 samples: the number of samples,
# the posterior probability that a sample in that group is a false positive,
# and the relative risk against the reference group. Both quantities come from
# the same posterior draws as the figure, so the table and Fig2.pdf can never
# disagree. Reference rows are kept here even where the figure omits them.

fmt_ci <- function(med, lo, hi, digits = 2) {
  ifelse(is.na(lo),
         sprintf(paste0("%.", digits, "f"), med),
         sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                 med, lo, hi))
}

tab_s_fig2 <- data.frame(
  covariate      = as.character(d_rr_full$covar),
  group          = as.character(d_rr_full$var),
  reference      = d_rr_full$ref,
  n              = d_rr_full$n,
  p_false_pos    = fmt_ci(d_rr_full$fp, d_rr_full$fp_lo, d_rr_full$fp_hi),
  rr_false_pos   = ifelse(as.character(d_rr_full$var) == d_rr_full$ref,
                          "1 (ref)",
                          fmt_ci(d_rr_full$rr, d_rr_full$rr_lo, d_rr_full$rr_hi)),
  stringsAsFactors = FALSE
)

# Numeric companion, for anyone who wants to recompute or re-plot.
tab_s_fig2_num <- d_rr_full[, c("covar", "var", "ref", "n",
                                "fp", "fp_lo", "fp_hi",
                                "rr", "rr_lo", "rr_hi")]
names(tab_s_fig2_num)[1:2] <- c("covariate", "group")

print(tab_s_fig2, row.names = FALSE)

write.csv(tab_s_fig2,     file.path(tab_dir, "TableS_Fig2_RR.csv"),     row.names = FALSE)
write.csv(tab_s_fig2_num, file.path(tab_dir, "TableS_Fig2_RR_raw.csv"), row.names = FALSE)

# Markdown version to paste straight into the supplement. The covariate name is
# printed once per block so the table reads the way the figure's facets look.
md_cov  <- ifelse(duplicated(tab_s_fig2$covariate), "", tab_s_fig2$covariate)
md_rows <- sprintf("| %s | %s | %d | %s | %s |",
                   md_cov, tab_s_fig2$group, tab_s_fig2$n,
                   tab_s_fig2$p_false_pos, tab_s_fig2$rr_false_pos)
writeLines(c(
  paste0("**Table S. Probability of a false positive diagnosis and relative ",
         "risk by covariate group, among samples with Ct<40 (Figure 2).** ",
         "Posterior medians with 95% credible intervals."),
  "",
  "| Covariate | Group | N | P(false positive) | RR vs reference |",
  "|---|---|---:|---|---|",
  md_rows
), file.path(tab_dir, "TableS_Fig2_RR.md"))



# =============================================================================
# Figure 3: Ct trajectory over time & serological results by model prob infection
# =============================================================================

# 3A: spline-fitted Ct trajectory; points coloured by P(true infection)
fig3a <- ggplot() +
  geom_point(data = pii, aes(t, ct, col = pi),
             alpha = 0.6, size = 1.2) +
  scale_color_distiller(palette = "RdYlBu", direction = 1,
                        name = "P(true\ninfection)",
                        limits = c(0, 1)) +
  geom_ribbon(data = tsi, aes(t, ymin = pred_lo, ymax = pred_hi),
              fill = "grey30", alpha = 0.08) +
  geom_ribbon(data = tsi, aes(t, ymin = mu_lo, ymax = mu_hi),
              fill = "grey30", alpha = 0.25) +
  geom_line(data = tsi, aes(t, mu_mean),
            col = "grey20", linewidth = 1) +
  theme_mpox + geom_vline(aes(xintercept=6), linewidth=1.5, linetype="dashed", alpha=0.9, col="seagreen")+
  labs(x = "Days since symptom onset", y = "Ct value")+ theme(legend.position = "none")
fig3a


# 3B serology figure
fig3b <- ggplot(sero[sero$ct>0, ], aes(days_symp, log_mfi, group=id)) +
  geom_point(aes(col=prob_inf, shape=factor(statusV2)), alpha=0.7) +
  geom_line(aes(col=prob_inf), alpha=0.2) +
  facet_wrap(~analyte, nrow=2) + theme_mpox +
  xlab("Days since symptom onset") +
  scale_color_distiller(palette = "RdYlBu", direction = 1,
                        name = "P(true\ninfection)",
                        limits = c(0, 1)) +
  scale_shape_manual(values = c(16, 17),
                     name = "Serostatus",
                     labels = c("Seronegative", "Seropositive")) +
  ylab("log(MFI)")+
  theme(legend.position="right")+ scale_x_continuous(trans="sqrt", breaks=c(0,10,100,300))
fig3b


# output
fig3 <- plot_grid(fig3a, fig3b, labels = c("A", "B"), rel_widths = c(1, 1.6))
ggsave(file.path(fig_dir, "Fig3.pdf"), fig3,
       width = 11, height = 5.5, units = "in", dpi = 300)


# =============================================================================
# Figure 4: ROC curve + case counts by Ct threshold
# =============================================================================

roc  <- ROC$rocq

# 4A: ROC curve
fig4a <- ggplot(roc[roc$cutoff>29,], aes(fpr, tpr, label = cutoff)) +
  geom_line(data=ROC$rocc, aes(fpr, tpr, group=iter), alpha=0.1, col="grey80")+
  geom_linerange(aes(ymin = tprL, ymax = tprU), linewidth = 0.5, col="blue") +
  geom_linerange(aes(xmin = fprL, xmax = fprU), linewidth = 0.5, col="blue") +
  geom_point(size = 2, col="blue") +
  geom_text(data = roc[roc$cutoff > 33, ],
            vjust = 3, hjust = - 0.5, col = "tomato3", size = 4.5) +
  theme_mpox + ylim(0.8,NA)+
  labs(x = "False positive rate", y = "Sensitivity")
fig4a

# 4B: TP / FP / FN case counts by threshold
df4b <- data.frame(
  cutoff = seq(30, 40),
  tp =  roc$tp[match(seq(30, 40), roc$cutoff)],
  fp =  roc$fp[match(seq(30, 40), roc$cutoff)],
  fn = -roc$fn[match(seq(30, 40), roc$cutoff)]
)
dff <- tidyr::pivot_longer(df4b, cols = c(tp, fp, fn), names_to = "group", values_to = "cases")
dff$group <- factor(dff$group, levels = c("fn", "fp", "tp"),
                    labels = c("False negative", "False positive",
                               "True positive"))
y_max <- max(dff$cases[dff$cases > 0])
y_min <- min(dff$cases[dff$cases < 0])

fig4b <- ggplot(dff, aes(factor(cutoff), cases, group = group, fill = group)) +
  geom_col(alpha = 0.8, col = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_fill_manual(values = c("indianred", "purple", "skyblue")) +
  scale_x_discrete(limits = rev) +
  coord_cartesian(clip = "off") +
  theme_mpox +
  theme(
    axis.title.y    = element_blank(),
    legend.title    = element_blank(),
    plot.margin     = margin(5.5, 5.5, 5.5, 40)
  ) +
  xlab("Ct threshold") +
  annotate("text", x = -Inf, y = 0.65 * y_max,
           label = "Detected cases",
           angle = 90, hjust = 0.5, vjust = -3.5,
           size = 3.8, inherit.aes = FALSE) +
  annotate("text", x = -Inf, y = 0.65 * y_min,
           label = "Undetected cases",
           angle = 90, hjust = 0.5, vjust = -3.5,
           size = 3.8, inherit.aes = FALSE)



# ===== Output =====
fig4 <- plot_grid(fig4a, fig4b, labels = c("A", "B"), rel_widths = c(0.7, 1))
ggsave(file.path(fig_dir, "Fig4.pdf"), fig4,
       width = 10, height = 5, units = "in", dpi = 300)

ROC$rocq$fnr  <- 1 - ROC$rocq$tpr
ROC$rocq$fnrL <- 1 - ROC$rocq$tprU
ROC$rocq$fnrU <- 1 - ROC$rocq$tprL


# =============================================================================
# Supp Figure: Empirical CDF of Ct values by site
# =============================================================================
fig_ctcdf <- ggplot(cts, aes(ct, colour = site)) +
  stat_ecdf(linewidth = 1) +
  theme_mpox +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Ct value", y = "Cumulative proportion < x", colour = "Site")

fig_ctcdf
ggsave(file.path(fig_dir, "Fig_CDF_Ct.pdf"), fig_ctcdf,
       width = 6, height = 4, units = "in", dpi = 300)


# =============================================================================
# Supp Figure: Site-specific infection probabilities
# =============================================================================
colnames(pars)[3:5] <- c("med","ciL","ciU")
site_ests <- ggplot(pars[pars$pars=="probI", ], aes(site, med))+ geom_point()+
  geom_linerange(aes(ymin=ciL, ymax=ciU))+ theme_mpox+ ylim(0,NA)+
  ylab("Prob(true infection | Ct<40)")+ theme(axis.title.x = element_blank())
site_ests


ggsave(file.path(fig_dir, "Fig_SiteProbI.pdf"), site_ests,
       width = 6, height = 4, units = "in", dpi = 300)


# =============================================================================
# Supp Figure - covariate effects on expected Ct values
# =============================================================================
baseline_ct <- tsi$mu_mean[which.min(abs(tsi$t - median(df$Tsymp_hosp)))]
ct_cri <- tsi[which.min(abs(tsi$t - median(df$Tsymp_hosp))),]


alpha_p <- alpha[alpha$group != "missing", ]
alpha_p$var   <- factor(alpha_p$var,   levels = covar_levels)
alpha_p$group <- factor(alpha_p$group,
                        levels = intersect(group_levels, unique(alpha_p$group)))

alpha_p$muL <- ct_cri$mu_lo
alpha_p$muU <- ct_cri$mu_hi

ct_covar <- ggplot(alpha_p, aes(group, exp(med) * baseline_ct)) +
  geom_hline(yintercept = baseline_ct,
             linetype = "dashed", col = "grey40", linewidth = 0.5) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=muL, ymax=muU), alpha=0.2, fill="grey")+
  geom_linerange(aes(ymin = exp(ciL) * baseline_ct,
                     ymax = exp(ciU) * baseline_ct),
                 col = "tomato", linewidth = 0.5) +
  geom_point(col = "tomato", size = 2) +
  facet_wrap(~var, scales = "free_x", nrow=2,
             labeller = label_wrap_gen(12)) +
  theme_mpox +
  theme(axis.text.x  = element_text(angle = 55, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Predicted Ct value")
ct_covar

# output SI fig
ggsave(file.path(fig_dir, "CTcovar_SIfig.pdf"), ct_covar,
       width = 10, height = 7, units = "in", dpi = 300)


# =============================================================================
# Supp Figure - ROC by platform
# =============================================================================
# read in data & model fits
stan_data  <- readRDS(file.path(res_dir, "Data.RDS"))
fit  <- readRDS(file.path(res_dir, "Fit.RDS"))

source(file.path(repo_dir, "Utils.R"))

# redo ROC analysis by platform
ROC_genexpert <- get_roc(fit, stan_data, cutoffs = seq(20, 40),
               negnoct = nrow(noamp[noamp$site %in% c("Uvira","Kamituga"), ]), cts_df = cts, sites = c("Uvira","Kamituga"))

ROC_radi <- get_roc(fit, stan_data, cutoffs = seq(20, 40),
                         negnoct = nrow(noamp[noamp$site %in% c("Goma","Kinshasa"), ]), cts_df = cts, sites = c("Goma","Kinshasa"))

# youdens J
ROC_genexpert$rocq$youden <- ROC_genexpert$rocq$tpr + (1 - ROC_genexpert$rocq$fpr) - 1
ROC_radi$rocq$youden <- ROC_radi$rocq$tpr + (1 - ROC_radi$rocq$fpr) - 1


# plot
rocplot1 <- ggplot(ROC_genexpert$rocq[ROC$rocq$cutoff > 29, ], aes(fpr, tpr, label = cutoff)) +
  theme_bw() +
  geom_line(data = ROC_genexpert$rocc, aes(fpr, tpr, group = iter), alpha = 0.1, col = "grey80") +
  xlim(0, NA) + labs(subtitle="Kamituga & Uvira (GeneXpert)")+
  geom_linerange(aes(ymin = tprL, ymax = tprU), col = "blue") +
  geom_point(col = "blue") +
  geom_linerange(aes(xmin = fprL, xmax = fprU), col = "blue") +
  geom_text(data=ROC_genexpert$rocq[ROC$rocq$cutoff > 33, ], vjust = 1.5, hjust = 1.3, col = "indianred2") +
  theme(text = element_text(size = 14), legend.position = "none") +
  xlab("False positive rate") + ylab("Sensitivity")+ ylim(0.8, NA)
rocplot1

rocplot2 <- ggplot(ROC_radi$rocq[ROC$rocq$cutoff > 29, ], aes(fpr, tpr, label = cutoff)) +
  theme_bw() +
  geom_line(data = ROC_radi$rocc, aes(fpr, tpr, group = iter), alpha = 0.1, col = "grey80") +
  xlim(0, NA) + labs(subtitle="Goma & Kinshasa (RADI)")+
  geom_linerange(aes(ymin = tprL, ymax = tprU), col = "blue") +
  geom_point(col = "blue") +
  geom_linerange(aes(xmin = fprL, xmax = fprU), col = "blue") +
  geom_text(data=ROC_radi$rocq[ROC$rocq$cutoff > 33, ], vjust = 1.5, hjust = 1.3, col = "indianred2") +
  theme(text = element_text(size = 14), legend.position = "none") +
  xlab("False positive rate") + ylab("Sensitivity")+ ylim(0.8, NA)
rocplot2

ROC_genexpert$rocq$assay <- "GeneXpert"
ROC_radi$rocq$assay <- "RADI"
roc_comp <- rbind(ROC_genexpert$rocq, ROC_radi$rocq)

tprassay <- ggplot(roc_comp[roc_comp$cutoff>29, ], aes(cutoff, tpr, group=assay, col=assay))+
  geom_point(position=position_dodge(width=0.5))+ ylab("Sensitivity")+ xlab("Ct threshold")+
  geom_linerange(aes(ymin=tprL, ymax=tprU), position=position_dodge(width=0.5))+ labs(col="")+
  theme_bw()+ theme(text=element_text(size=14), legend.position = c(0.7,0.3))+ ylim(0.75,NA)
fprassay <- ggplot(roc_comp[roc_comp$cutoff>29, ], aes(cutoff, fpr, group=assay, col=assay))+
  geom_point(position=position_dodge(width=0.5))+ ylab("False positive rate")+ xlab("Ct threshold")+
  geom_linerange(aes(ymin=fprL, ymax=fprU), position=position_dodge(width=0.5))+
  theme_bw()+ theme(text=element_text(size=14), legend.position = "none")

plot_grid(tprassay, fprassay)
rocbyassay1 <- plot_grid(rocplot1, rocplot2)
rocbyassay2 <- plot_grid(tprassay, fprassay)
rocbyassay <- plot_grid(rocbyassay1, rocbyassay2, ncol=1, labels=c("A","B"))

# output SI fig
ggsave(file.path(fig_dir, "ROC_byassay.pdf"), rocbyassay,
       width = 8, height = 6, units = "in", dpi = 300)


