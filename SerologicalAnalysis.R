#-------- Script to fit mixture models to serological data & align with Ct results -----------------#

# libraries
library(stringr)
library(ggplot2)
library(posterior)
library(bayesplot)
library(cmdstanr)
library(cowplot)
library(GGally)

# ===== Paths =====
repo_dir <- getwd()
data_dir <- file.path(repo_dir, "data")
stan_dir <- file.path(repo_dir, "stan")
out_dir  <- file.path(repo_dir, "results/serology")
fig_dir  <- file.path(repo_dir, "results/figures")
dir.create(out_dir, showWarnings = FALSE)

# ---- read in serology data ------------------------------------------

dfs <- read.csv(file.path(data_dir, "Serology.csv"))

# ---- QC: drop failed assay wells ------------------------------------------
# Whole samples (all 8 analytes at once) that break sharply from the
# individual's own trajectory and recover at the following visit. Each is a
# well-level assay failure rather than a biological signal:
#      1 S6 - every analyte at the global minimum of its distribution, 100-1000x
#             below the flanking visits; blank-well signature.
#     74 S6 - same signature (e.g. MPXV A27 875 -> 19 -> 1011).
#     82 S6 - mirror image: flat baseline for five visits, all 8 analytes jump
#             20-140x, then return at S7. The profile is a near-clone of
#             participant 1's S7 sample (cor of logs 0.99), so 1 S6 and 82 S6 look
#             like a single aliquot mix-up on the day-270 run.
# All three are post-baseline visits, so serostatus (defined from S1/S2) and
# the analytic N are unchanged; the drops affect the mixture fits and figures.
# these will be retested but not on timescale for this paper
qc_drop <- data.frame(id    = c(1, 74, 82),
                      visit = c("S6",  "S6",  "S6"),
                      stringsAsFactors = FALSE)
dfs <- dfs[!(paste(dfs$id, dfs$visit) %in% paste(qc_drop$id, qc_drop$visit)), ]

# ---- build Stan input list per antigen -----------------------------------
ags <- unique(dfs$analyte)

stan_data <- setNames(lapply(ags, function(a) {
  dfa   <- dfs[dfs$analyte == a, ]
  y     <- log(dfa$mfi)
  y_fit <- seq(min(y) - 0.5, max(y) + 0.5, 0.2)
  list(N = nrow(dfa), y = y, y_fit = y_fit, predL = length(y_fit))
}), ags)

source(file.path(repo_dir, "Utils.R"))

# ---- compile Stan model --------------------------------------------------
check_cmdstan_toolchain()
# set_cmdstan_path("path/to/cmdstan-X.X.X")
mod <- cmdstan_model(file.path(stan_dir, "SeroMixture.stan"), pedantic = FALSE)
color_scheme_set("mix-blue-red")

# ---- fit mixture model to each antigen independently ---------------------
traceplots    <- list()
par_ests      <- list()
mixture_fits  <- list()
serostatus    <- list()

# Fixed seeds so the eight mixture fits are reproducible run to run. Both
# sources of randomness have to be pinned: Stan's sampler via `seed`, and the
# R-side starting values, since make_sero_inits() draws from runif(). Without
# this the per-antigen serostatus calls wobble by one or two between runs
# (the universal serostatus has been stable, but do not rely on that).
sero_seed <- 20260821

for (i in seq_along(ags)) {
  a <- ags[i]
  set.seed(sero_seed + i)          # governs make_sero_inits() below
  fit <- mod$sample(
    data            = stan_data[[a]],
    chains          = 3,
    parallel_chains = 3,
    iter_sampling   = 2000,
    iter_warmup     = 1000,
    refresh         = 100,
    seed            = sero_seed + i,
    init            = make_sero_inits(3)
  )
  
  draws             <- fit$draws(format = "df")
  traceplots[[a]]   <- mcmc_trace(draws, regex_pars = c("sero", "mu", "sd", "lp__")) + labs(subtitle = a)
  par_ests[[a]]     <- extract_mixture_pars(draws)
  mixture_fits[[a]] <- plot_mixture_fits(fit, stan_data[[a]], a)
  serostatus[[a]]   <- extract_serostatus(fit, dfs[dfs$analyte == a, ])
}

#--- check convergence & fits

# traceplots
plot_grid(plotlist = traceplots) 

# plot model fits
shared_legend      <- get_legend(mixture_fits[[1]])
mixture_fits_clean <- lapply(mixture_fits, function(p) {
  p + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
})

inner_grid <- plot_grid(plotlist = mixture_fits_clean, nrow = 2)
with_x     <- plot_grid(inner_grid,
                        ggdraw() + draw_label("log(MFI)", size = 14),
                        ncol = 1, rel_heights = c(1, 0.04))
with_xy    <- plot_grid(ggdraw() + draw_label("Density", angle = 90, size = 14),
                        with_x,
                        ncol = 2, rel_widths = c(0.04, 1))
plotfits <- plot_grid(with_xy, shared_legend, rel_widths = c(1, 0.12))
plotfits


# ---- define serostatus ----------------------------------------------------

# compile results
serostatus_df <- do.call(rbind, Map(function(df, a) { df$analyte <- a; df }, serostatus, names(serostatus)))
rownames(serostatus_df) <- NULL
serostatus_df$visitT <- factor(serostatus_df$visit, labels=c("0","14","30","90","180","270","360"))


# per-antigen: seropositive if (a) component prob > 0.5 at visit 0 or 14,
#              OR (b) four-fold rise in MFI between visit 0 and visit 14
# Minimum MFI guard: the prob>0.5 criterion only fires if the individual's
# log(MFI) is at or above the estimated seronegative component mean (mu0).
# This prevents the broader seropositive component from "collecting" outlier
# observations that fall below the entire seronegative distribution.
mu0_by_analyte <- sapply(par_ests, function(p) p$med[p$par == "mu0"])

v12 <- serostatus_df[serostatus_df$visitT %in% c("0", "14"), ]
v12$seropos_adj <- v12$seropos & v12$log_mfi >= mu0_by_analyte[v12$analyte]
v12_ag_status <- aggregate(seropos_adj ~ id + analyte, data = v12, FUN = any)
names(v12_ag_status)[3] <- "statusV2_ag"

# four-fold rise: log_mfi_v14 - log_mfi_v0 >= log(4)
v0   <- v12[v12$visitT == "0",  c("id", "analyte", "log_mfi")]
v14  <- v12[v12$visitT == "14", c("id", "analyte", "log_mfi")]
rise <- merge(v0, v14, by = c("id", "analyte"), suffixes = c("_v0", "_v14"))
rise$fourfold <- rise$log_mfi_v14 - rise$log_mfi_v0 >= log(4)
v12_ag_status <- merge(v12_ag_status, rise[, c("id", "analyte", "fourfold")],
                       by = c("id", "analyte"), all.x = TRUE)
v12_ag_status$fourfold[is.na(v12_ag_status$fourfold)] <- FALSE
v12_ag_status$pos_criterion <- ifelse(
  v12_ag_status$statusV2_ag & v12_ag_status$fourfold, "Both",
  ifelse(v12_ag_status$fourfold,                       "Four-fold rise only",
  ifelse(v12_ag_status$statusV2_ag,                    "Prob > 0.5 only", "Negative"))
)
v12_ag_status$statusV2_ag <- v12_ag_status$statusV2_ag | v12_ag_status$fourfold

# save criterion lookup separately for plotting; drop working columns
criterion_lookup <- v12_ag_status[, c("id", "analyte", "pos_criterion")]
v12_ag_status    <- v12_ag_status[, c("id", "analyte", "statusV2_ag")]

serostatus_df <- merge(serostatus_df, v12_ag_status, by = c("id", "analyte"), all.x = TRUE)

# universal: seropositive if a majority of the eight antigens responded.
# The per-participant count of responding antigens is sharply bimodal -- most
# respond on none and the next largest group on all eight, with only a thin
# scatter in between, almost all of those on one or two antigens -- so the
# majority cut falls in a natural gap rather than slicing through a continuum.
# That intermediate group is not a partial serological response: its median Ct
# and its median model-inferred P(true infection) are indistinguishable from
# the fully seronegative group (roughly Ct 39 and P=0), whereas the
# majority-positive group sits near Ct 18 and P=1. Requiring only one antigen
# out of eight gives the assay eight chances to cross prob 0.5 and admits that
# noise. ManuscriptNumbers.R prints the current counts for this comparison.
# (alternative: any antigen positive — any(x, na.rm = TRUE))
ag_counts  <- unique(serostatus_df[, c("id", "analyte", "statusV2_ag")])
v12_status <- aggregate(statusV2_ag ~ id, data = ag_counts, FUN = function(x) mean(x, na.rm = TRUE) > 0.5)
names(v12_status)[2] <- "statusV2"
serostatus_df <- merge(serostatus_df, v12_status, by = "id", all.x = TRUE)

# compare per-antigen vs universal serostatus classification
status_compare <- unique(serostatus_df[, c("id", "analyte", "statusV2", "statusV2_ag")])
table(universal = status_compare$statusV2, per_antigen = status_compare$statusV2_ag)

# save results
serostatus_df$statusV2 <- factor(serostatus_df$statusV2, labels=c("Negative","Positive"))
serostatus_df$statusV2_ag <- factor(serostatus_df$statusV2_ag, labels=c("Negative","Positive"))
serodf <- merge(serostatus_df, dfs[,c("id","visit","analyte","days_symp","ct","prob_inf")], by = c("id", "visit", "analyte"))
keep_ids <- with(serodf, names(which(
  tapply(visitT, id, function(v) any(c("0", "14") %in% v))
)))
serodf <- serodf[serodf$id %in% keep_ids, ]
serodf <- merge(serodf, criterion_lookup, by = c("id", "analyte"), all.x = TRUE)

# output
write.csv(serodf, file.path(out_dir, "SeroResults.csv"), row.names = FALSE)


#----- plot results -------------------------------------------------------

# plot trajectories by universal serostatus
fig_seroclass <- ggplot(serodf, aes(days_symp, log_mfi))+ geom_point(aes(col=statusV2))+ 
  geom_line(aes(group=id), alpha=0.3, col="grey")+
  facet_wrap(~analyte, ncol=4)+ theme_bw()+ labs(col="Overall Response")+ scale_color_manual(values=c("skyblue","orange"))+
  ylab("log(MFI)")+ xlab("Days since symptom onset")+ theme(text=element_text(size=14))
fig_seroclass

# plot by any antigen serostatus
fig_seroclass <- ggplot(serodf, aes(days_symp, log_mfi))+ geom_point(aes(col=statusV2_ag))+
  geom_line(aes(group=id), alpha=0.3, col="grey")+
  facet_wrap(~analyte, ncol=4)+ theme_bw()+ labs(col="Overall Response")+ scale_color_manual(values=c("skyblue","orange"))+
  ylab("log(MFI)")+ xlab("Days since symptom onset")+ theme(text=element_text(size=14))
fig_seroclass

# plot trajectories by how seropositivity was determined (four-fold rise vs prob>0.5)
crit_cols <- c("Negative" = "grey70", "Prob > 0.5 only" = "orange",
               "Four-fold rise only" = "steelblue", "Both" = "purple")
serodf$pos_criterion <- factor(serodf$pos_criterion,
  levels = c("Negative", "Prob > 0.5 only", "Four-fold rise only", "Both"))
fig_criterion <- ggplot(serodf, aes(days_symp, log_mfi, group = id)) +
  geom_line(aes(col = pos_criterion), alpha = 0.4) +
  geom_point(aes(col = pos_criterion), alpha = 0.7) +
  scale_color_manual(values = crit_cols, name = "Classification") +
  facet_wrap(~analyte, ncol = 4) +
  theme_bw() + theme(text = element_text(size = 14)) +
  ylab("log(MFI)") + xlab("Days since symptom onset")
fig_criterion

# summary of responses
serosumm <- serodf[!duplicated(serodf$id), ] # 1 entry per person for summaries
table(serosumm$statusV2)
table(serosumm$statusV2[serosumm$ct>0])


# summary of prob true infection from Ct model by serostatus
summary(serosumm$prob_inf[serosumm$statusV2=="Negative" & serosumm$ct>0])
summary(serosumm$prob_inf[serosumm$statusV2=="Positive"  & serosumm$ct>0])
ggplot(serosumm, aes(statusV2, prob_inf))+ geom_boxplot()+ geom_jitter(width=0.1)+
  theme_bw()+ xlab("Serological status")+ ylab("P(true inf | Ct<40) from Ct model")

# summary of Ct values by serostatus
summary(serosumm$ct[serosumm$statusV2=="Negative" & serosumm$ct>0])
summary(serosumm$ct[serosumm$statusV2=="Positive"  & serosumm$ct>0])

# serostatus among inferred false positive - single antigen
table(serosumm$statusV2_ag[serosumm$ct<40 & serosumm$prob_inf<0.5])
table(serosumm$statusV2_ag[serosumm$ct<40 & serosumm$prob_inf>=0.5])
table(serosumm$statusV2_ag[serosumm$ct==0])

# concordance
# need to exclude people with no cts
concordant <- sum(
  (serosumm$ctfp == 0 & serosumm$statusV2_ag == "Positive") |
    (serosumm$ctfp == 1 & serosumm$statusV2_ag == "Negative")
)
percent_agreement <- concordant / nrow(serosumm)

library(irr)
serosumm$serostat <- as.numeric(as.factor(serosumm$statusV2_ag))-1
kappa2(serosumm[, c("ct", "serostat")])

# plot serological trajectories by Ct values
fig_seroct <- ggplot(serodf[serodf$ct>0, ], aes(days_symp, log_mfi))+ geom_point(aes(col=ct))+ 
  geom_line(aes(group=id, col=ct), alpha=0.3)+ scale_color_viridis_c()+
  facet_wrap(~analyte, ncol=4)+ theme_bw()+ labs(col="Ct value")+ 
  ylab("log(MFI)")+ xlab("Days since symptom onset")+ theme(text=element_text(size=14))+
  ggnewscale::new_scale_colour() +
  geom_point(data=serodf[serodf$ct==0, ], aes(colour = "Negative (no Ct)"), alpha=0.6) +
  geom_line(data=serodf[serodf$ct==0, ], aes(group=id, colour = "Negative (no Ct)"), alpha=0.6) +
  scale_colour_manual(values = c("Negative (no Ct)" = "grey"), name = NULL)
fig_seroct


# plot serological trajectories by Ct model inferences
serodf$group <- "True positive"
serodf$group[serodf$prob_inf<0.5] <- "False positive"
serodf$group[serodf$ct==0 | serodf$ct>=40] <- "Ct>=40 or no amplification"

fig_seroctmod <- ggplot(serodf, aes(days_symp, log_mfi))+ geom_point(aes(col=group))+ 
  geom_line(aes(group=id, col=group), alpha=0.3)+ 
  facet_wrap(~analyte, ncol=4)+ theme_bw()+ labs(col="Ct model inference")+ 
  ylab("log(MFI)")+ xlab("Days since symptom onset")+ theme(text=element_text(size=14))
fig_seroctmod  


# distribution of cts by serological status
fig_seroctclass <- ggplot(serodf[!duplicated(serodf$id) & serodf$ct>0, ], aes(statusV2_ag, ct, fill=statusV2_ag))+ theme_bw()+
  geom_violin()+ geom_jitter(width=0.02)+ theme(text=element_text(size=14), legend.position="none")+
  ylab("Ct value")+ xlab("Overall serological response")+ scale_fill_manual(values=c("skyblue","orange"))
fig_seroctclass

#--------- output supp figures
ggsave(file.path(fig_dir, "SeroFits.pdf"),           plotfits, width=11, height=5, units="in", dpi=300)
ggsave(file.path(fig_dir, "SeroClassification.pdf"), fig_seroclass, width=11, height=6, units="in", dpi=300)
ggsave(file.path(fig_dir, "SeroCt.pdf"),             fig_seroct,    width=11, height=6, units="in", dpi=300)
ggsave(file.path(fig_dir, "SeroCtClass.pdf"),        fig_seroctclass, width=5, height=5, units="in", dpi=300)
ggsave(file.path(fig_dir, "SeroCtModelGroup.pdf"),   fig_seroctmod, width=11, height=6, units="in", dpi=300)


