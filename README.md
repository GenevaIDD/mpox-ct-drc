# Mpox Ct value analysis

Bayesian latent class model to estimate the probability that a positive MPXV qPCR result from a suspected mpox case represents a true infection, using Ct values from four sites in the Democratic Republic of Congo. This repository reproduces the main analyses in O'Driscoll et al., 2026, *"The risk of mpox false positive results in high transmission settings: evidence from a multi-site observational study in DR Congo"*, Lancet Infectious Diseases, [doi:10.1016/S1473-3099(26)00411-1](https://doi.org/10.1016/S1473-3099(26)00411-1).

## Model

The main model (`stan/LatentClassCtModel.stan`) fits a two-component latent class model to log-transformed Ct values:

-   **Infection component** — right-truncated Normal(μ, σ_I) where μ follows a B-spline trajectory over days since symptom onset plus covariate offsets, truncated at a specified Ct threshold
-   **Environmental component** — Normal(μ_E, σ_E) with parameters anchored to environmental qPCR sample results

The mixture weight π (proportion of true infections) is estimated per site. Covariate effects on Ct (age, sex, sample type, HIV status, vaccination, symptom profiles) are modelled with hierarchical shrinkage.

## Requirements

**R packages:** `cmdstanr`, `splines`, `posterior`, `loo`, `bayesplot`, `ggplot2`, `cowplot`, `dplyr`, `ggrepel`, `GGally`, `matrixStats`, `stringr`, `ggh4x`, `ggridges`, `irr`, `rnaturalearth`, `rnaturalearthdata`, `sf`

**Stan:** Install via `cmdstanr::install_cmdstan()`.

## Repository structure

```         
FitMainModel.R          # Fit the main latent class model and save results
SerologicalAnalysis.R   # Fit serology mixture models and classify serostatus
MakeFigures.R           # Reproduce manuscript figures from saved results
Descriptives.R          # Compute descriptive statistics reported in the Results
ManuscriptNumbers.R     # Recompute every number in the paper and check it against the published value
LoadCaseCts.R           # Reader for data/CaseCts.csv, sourced by every script that uses it
Utils.R                 # Helper functions sourced by all fitting scripts
stan/
  LatentClassCtModel.stan   # Main latent class model
  NullCtModel.stan          # Null model, no environmental component
  SeroMixture.stan          # Two-component mixture model for serological data
data/
  CaseCts.csv               # Case Ct values and covariates
  EnvironmentalSamplingCts.csv  # Environmental qPCR Ct values
  Serology.csv              # Serological data
results/
  main/                     # Main model outputs (fit objects and summary CSVs)
  serology/                 # Serostatus classifications (SeroResults.csv)
  figures/                  # All output figures
  tables/                   # Supplementary tables (regenerated, not tracked)
```

All three input files are plain CSVs, so they can be inspected without R. Read the case data with `read_case_cts()` from `LoadCaseCts.R` rather than `read.csv()` directly: it restores the natural level order of `age_group` (`0-4`, `5-9`, `10-14`, ...), which determines both the integer covariate coding sent to Stan and the labels applied to the model output. A bare `read.csv()` leaves `age_group` as a character vector that later sorts alphabetically (`0-4`, `10-14`, `15-19`, ...), which silently mislabels the age effects.

## Usage

All scripts should be run from the repository root. The recommended order is:

### 1. Descriptive statistics

``` r
source("Descriptives.R")
```

Prints all participant-level and Ct distribution summary statistics reported in the Results section. Runs in seconds; no model fitting required.

Two numbers from the manuscript cannot be **recomputed** from the public data:

-   **Median age (IQR)** — the Table 2 row is **hard-coded** in `Descriptives.R` from the original patient-level dataset. Only `age_group` is shared in `CaseCts.csv`, so numeric age is not available here; the printed values are correct but will not change if you modify the data.
-   **N=243 paired lesion/oropharyngeal swabs** — de-duplication was applied upstream before the data were shared, so this count is not recoverable and is not printed.

`ManuscriptNumbers.R` (below) lists every such gap across the whole paper, not just this section.

### 2. Fit the main model

``` r
source("FitMainModel.R")
```

Fits the latent class Stan model (~10–30 min depending on hardware) and saves all results to `results/main/`. Requires CmdStan — install via `cmdstanr::install_cmdstan()`.

### 3. Fit the serological mixture models

``` r
source("SerologicalAnalysis.R")
```

Fits a two-component Normal mixture to each serological antigen independently and classifies participants as seropositive or seronegative. Saves results to `results/serology/`.

### 4. Reproduce figures

``` r
source("MakeFigures.R")
```

Reads saved CSVs and the full fit object to produce all manuscript figures. **Requires `results/main/Fit.RDS`**, which is not committed to the repository due to file size — run `FitMainModel.R` first to generate it.

### 5. Check every number in the paper

``` r
source("ManuscriptNumbers.R")
```

Recomputes every figure quoted in the paper — abstract, both tables, all four figures, and the in-text estimates — from `data/` and the saved outputs in `results/main`, and prints each one next to the value in the published proofs with a verdict of `MATCH`, `DIFFERS`, or `NO DATA`. The audit is also written to `results/ManuscriptNumbers.csv`.

The published values are those in the fourth set of proofs (26TLID0924_Azman, saved 09-Sept-2026). A closing section answers the `[A: ...]` author queries in that proof.

Base R only, and no refitting: it reads the saved CSVs and `ROC.RDS`. Two optional extras:

-   Set `MPOX_RESULTS_DIR` to audit a different output directory, e.g. `MPOX_RESULTS_DIR=results/main_refit Rscript ManuscriptNumbers.R`.
-   Credible intervals for the figure 2 relative risks are ratios of two group probabilities, so they have to be formed draw by draw and need `results/main/Fit.RDS`. When it is present the script computes them once and caches them as `RR_covars.csv` alongside it; without it the relative risks are printed as point estimates only.

The script recovers each covariate group label from `Data.RDS` rather than trusting the labels in `pI_covars.csv`, so results produced before the August 2026 covariate-coding fix (which left `age_group` labelled in alphabetical order) are still reported against the right groups. It warns when it finds such a file.

## Outputs

`FitMainModel.R` saves the following to `results/main/`:

| File | Contents |
|-----------------------|-------------------------------------------------|
| `Fit.RDS` | Full `cmdstanr` fit object |
| `Data.RDS` | Stan data list passed to the model |
| `Params.csv` | Posterior summaries for key parameters |
| `pii.csv` | Per-individual P(true infection) |
| `pI_covars.csv` | P(true infection \| Ct \< threshold) by covariate group |
| `covar_coef.csv` | Covariate effects on Ct |
| `TimeSinceSymps.csv` | Spline-fitted Ct trajectory over time |
| `FitPDFs.csv` | Population-marginal mixture component densities |
| `WAICLOO.csv` | LOO and WAIC estimates |
| `ROC.RDS` | ROC analysis across Ct thresholds |

`MakeFigures.R` saves the manuscript figures to `results/figures/` and, alongside figure 2, the supplementary table behind it to `results/tables/`:

| File | Contents |
|-----------------------|-------------------------------------------------|
| `TableS_Fig2_RR.csv` | Formatted figure 2 table: N, P(false positive), and RR by covariate group |
| `TableS_Fig2_RR_raw.csv` | Same rows with unrounded posterior medians and interval bounds |
| `TableS_Fig2_RR.md` | Markdown version to paste into the supplement |

The probabilities and relative risks in the table come from the same posterior draws as figure 2, so the two cannot drift apart. Reference groups are kept as rows (`RR = 1 (ref)`) even where the figure drops them, and groups with a missing covariate value are excluded, as in the figure.

`SerologicalAnalysis.R` saves `results/serology/SeroResults.csv` (per-visit serostatus and MFI values) and supplementary figures to `results/figures/`.

`ManuscriptNumbers.R` saves `results/ManuscriptNumbers.csv` (one row per checked number: section, item, published value, recomputed value, verdict, note).

## Participant identifiers

The `id` column in `data/Serology.csv` and `results/serology/SeroResults.csv` is a random study-specific identifier assigned for this release. It does not correspond to any identifier used in the field or laboratory, and it is not linked to the `id` column in `data/CaseCts.csv`.

## License

The code in this repository is licensed under the GNU General Public License v3.0 (see `LICENSE`). The data in `data/` and the model outputs in `results/` are licensed under the Creative Commons Attribution 4.0 International license (see `LICENSE-DATA`); please cite O'Driscoll et al., 2026 when reusing them.
