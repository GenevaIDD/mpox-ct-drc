# =============================================================================
# Case Ct data loader
# =============================================================================
# Single entry point for data/CaseCts.csv, sourced by every script that reads
# the case data. Base R only, no package dependencies.
#
# The file is a plain CSV so that it can be inspected without R, which means the
# categorical covariates arrive as character vectors. Only age_group has a
# meaningful order that is not the alphabetical one, and that order is
# load-bearing: FitMainModel.R derives the integer codes it sends to Stan with
# levels(factor(cts[[v]])), so a level order of 0-4, 10-14, 15-19, ... would
# silently shift both the covariate coding and the labels applied to the model
# output. read_case_cts() restores the intended order on every read.

AGE_GROUP_LEVELS <- c("0-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40+")

read_case_cts <- function(data_dir) {
  df <- read.csv(file.path(data_dir, "CaseCts.csv"), stringsAsFactors = FALSE)

  unknown <- setdiff(unique(df$age_group), AGE_GROUP_LEVELS)
  if (length(unknown))
    stop("age_group values not in AGE_GROUP_LEVELS: ", paste(unknown, collapse = ", "))
  df$age_group <- factor(df$age_group, levels = AGE_GROUP_LEVELS)

  df
}
