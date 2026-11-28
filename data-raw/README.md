# data-raw: simulated critical values and p-values

Reproduces `R/sysdata.rda`. Run everything from this folder; the simulations
need Rcpp. Nothing here is part of the built package (`.Rbuildignore`).

## Engle-Granger residual-based ADF statistic: `eg_rs`

1. `eg_simulate.R` (kernel `tau_kernel.cpp`) — the statistic under the null
   of no cointegration: `nvars` independent random walks, Step 1 regression
   of the first on the others with no deterministic terms (`n`), a constant
   (`c`) or a constant and trend (`ct`), then the ADF t ratio of the
   residuals without deterministic terms on the maximal sample, as
   `unitroottools::adf(trend = "n")`. Same design as unitroottools' τ
   simulation (D21): T = 20, ..., 1000; k = 0, ..., min(12, Schwert's kmax),
   plus k = 13, ..., 21 for T >= 150; 40 experiments x 10,000 replications
   per (nvars, T); quantiles at 117 probabilities. Seeds 20261126 + nvars
   (L'Ecuyer-CMRG streams over T).

   ```bash
   cd data-raw
   Rscript eg_simulate.R          # nvars = 2..6
   Rscript eg_simulate.R 7 8      # more nvars later, same design
   ```

   One job per (nvars, T), largest first, on all cores. Each finished job is
   saved in `eg_parts/` (not committed); if the run stops, the same command
   resumes it. At the end it writes `eg_quantiles_n<nvars>.csv.gz`.

   Cost: about 2 hours per nvars on 2 cores of the development machine
   (timed on a small run), so about 20 core-hours for nvars = 2..6.

2. `Rscript eg_response_surfaces.R` — `fit_surface()` (from `surface_fit.R`)
   for each nvars and case on T >= 25: a common asymptotic value plus a cubic
   in 1/T for each k. Writes `eg_rs.rds`, reports the fit, checks that the
   fitted quantiles are monotone, and compares the k = 0 critical values with
   MacKinnon (1996) (his tables as shipped in urca; urca is only needed
   here).
3. `Rscript make_sysdata.R` — writes `../R/sysdata.rda` (`eg_rs`, and `ecm_rs` once it exists).
4. `Rscript eg_validate.R` — size of the 1/5/10% critical values and of
   p <= level in fresh simulations at (T, k) not on the grid
   (T = 33, 60, 137, 290, 610; k up to 20), 20,000 replications per cell.

## ECM test statistic: `ecm_rs`

1. `ecm_simulate.R` (kernel `ecm_kernel.cpp`) — the t ratio of y_{t-1} in
   the unrestricted ECM, dy_t on [deterministics, y_{t-1}, x_{t-1}, dx_t,
   dy_{t-1..t-p}, dx_{t-1..t-p}], plain OLS standard error, under the null
   (`nvars` independent random walks). Cases `c` and `ct`; p = 0..4 lagged
   differences; T = 20, ..., 1000; 40 x 10,000 replications per (nvars, T);
   seeds 20261127 + nvars. Resumable through `ecm_parts/` (not committed).

   ```bash
   cd data-raw
   Rscript ecm_simulate.R          # nvars = 2..6
   ```

   Cost: timed on a small run, about 9 core-hours of the development
   machine for nvars = 2..6, so under an hour on 8-10 cores.
2. `Rscript ecm_response_surfaces.R` — surfaces (p in the role of k),
   monotonicity check, and the 1/5/10% values by p for comparison with
   Ericsson and MacKinnon (2002) and the chapter's Table 2.
3. `Rscript make_sysdata.R` — now stores `eg_rs` and, if present, `ecm_rs`.
4. `Rscript ecm_validate.R` — size in fresh simulations at T = 33, 60, 137,
   290, 610 and p up to 4.

## Shared code

`tau_kernel.cpp` and `surface_fit.R` are copies of the files in
unitroottools/data-raw, and `../R/calibration.R` of unitroottools'
`R/calibration.R`; keep them identical. For nvars = 1 the kernel gives the
plain ADF τ statistics, i.e. unitroottools' own tables.

## Phillips-Ouliaris Z_t and Z_alpha: `po_rs`

1. `Rscript po_simulate.R` (kernel `po_kernel.cpp`) — the two statistics
   under the null, with the Engle-Granger DGP and cases (`n`, `c`, `ct`),
   nvars = 2..6, T = 20..1000, 40 x 10,000 replications per (nvars, T);
   Bartlett kernel with Andrews' bandwidth, the defaults of `po_test()`. No
   lag order (k = 0). Seeds 20261128 + nvars. Writes `po_quantiles.csv.gz`.
   About 20 minutes on 2 cores (it was run in the development sandbox).
2. `Rscript po_response_surfaces.R` — Z_t: response surface as Engle-Granger;
   Z_alpha: tabulated by T (`table_calibration()`), because its quantiles are
   not a cubic in 1/T. Writes `po_rs.rds` and the log `po_surfaces.log`.
3. `Rscript make_sysdata.R`.
4. `Rscript po_validate.R` — size of the 5% tests in fresh simulations
   (`po_validate.log`).

## Gregory-Hansen tables

The asymptotic critical values are typed into `R/gh_test.R` from Gregory and
Hansen (1996a, b). `Rscript gh_check.R` compares them with a simulation at
T = 400 (`gh_check.log`; a few minutes on 2 cores).
