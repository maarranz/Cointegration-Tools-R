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
3. `Rscript make_sysdata.R` — writes `../R/sysdata.rda`.
4. `Rscript eg_validate.R` — size of the 1/5/10% critical values and of
   p <= level in fresh simulations at (T, k) not on the grid
   (T = 33, 60, 137, 290, 610; k up to 20), 20,000 replications per cell.

## Shared code

`tau_kernel.cpp` and `surface_fit.R` are copies of the files in
unitroottools/data-raw, and `../R/calibration.R` of unitroottools'
`R/calibration.R`; keep them identical. For nvars = 1 the kernel gives the
plain ADF τ statistics, i.e. unitroottools' own tables.
