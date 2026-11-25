# cointegrationtools — design notes (planning draft)

Planning draft, 25 November 2026. Repository `maarranz/Cointegration-Tools-R`,
package `cointegrationtools`. Like the `DESIGN.md` of adltools, this file is
excluded from the built package (`.Rbuildignore`) and records the reasons
behind the API. Items marked **[open]** are not decided yet.

## 1. Place in the ATSE suite

| Package | Role |
|:--|:--|
| arimatools | ARIMA estimation, diagnostics, outliers, dynamic forecasting |
| forecastingtools | Static forecasts, accuracy, evaluation and comparison |
| unitroottools | Unit-root and stationarity tests, deterministic-term strategies |
| adltools | Stationary ADL/ECM models: estimation, inference, Granger, diagnostics, stability, exogeneity, dynamic forecasts |
| **cointegrationtools** | Single-equation cointegration: cointegrating regression, residual-based and ECM-based tests, the cointegrated ECM, and (later) breaks in the cointegrating vector |

- **Dependency direction** (agreed in adltools' DESIGN.md §1):
  cointegrationtools → adltools, unitroottools; never the reverse.
- **Out of scope: Johansen / VECM.** The chapter treats it as reading
  material, and `urca::ca.jo()` covers it. The package handles any number of
  regressors **with one cointegrating relation** (money demand being the
  standard case).
- **In scope eventually, by adltools' own list** (adltools §1, "Cointegration
  stays out"): bounds tests, ECM-based cointegration tests, FM-OLS/DOLS,
  Gregory–Hansen. Bounds test (Pesaran–Shin–Smith) **[open]**: see §11.

## 2. Notation (as in Cointegration_T)

| Symbol | Meaning |
|:--|:--|
| y_t, x_{i,t} | dependent variable and the n − 1 regressors; n = number of I(1) variables in the relation |
| α, β (β_i) | intercept and cointegrating coefficients of the levels regression y_t = α + β′x_t + z_t |
| z_t, ẑ_t | equilibrium error and its estimate |
| δ₀, δ_i, γ_j | ECM intercept, lagged Δy and Δx coefficients |
| b | speed of adjustment (coefficient on ẑ_{t−1}) |
| θ, π | ECM test: coefficient on (y − x)_{t−1} (or y_{t−1}) and on x_{t−1} |
| τ* , τ̂* | Gregory–Hansen break fraction (later) |

Printed output uses these labels, as adltools uses α/β/k/b from ADL1–ADL4.

## 3. Functions by stage

### 3.1 Pre-test: integration order — `integration_order()`

A table, one row per series, from `unitroottools::ur_strategy()` (or `adf()`
with user-fixed deterministics), plus the same on first differences.
Answers Step 1 of the chapter's applied example ("is each series I(1)?") in
one call. Thin wrapper, no new testing code. Warns if any series looks I(0)
or I(2) (the CI(2,1) case is outside the tests' tables; see §10).

### 3.2 Cointegrating regression — `coint_reg()`

```r
coint_reg(y ~ x1 + x2, data,
          method = c("OLS", "FM", "D", "IM"),
          trend  = c("c", "ct", "n"),
          kernel = "ba", bandwidth = "and",            # FM, IM inference
          leads = NULL, lags = NULL, info_crit = "AIC") # DOLS
```

- **Engine: own code** (decided 25 Nov 2026; `R/lrv.R`, `R/estimators.R`).
  cointReg (CRAN 0.2.0, 2016) is in Suggests only, as the reference in the
  tests. Reasons, from reading its source:
  - `cointRegIM()` returns `varmat` without ω̂_{u·v}: used as a covariance
    it overstates the standard errors (≈50% in a test case).
  - `cointRegD()` compares AIC/BIC across lead/lag pairs fitted on
    different samples, with a penalty that ignores the deterministic terms;
    with 0 leads and 0 lags it drops Δx_t and returns static OLS.
  - p-values from t(T − k); coefficient names `deter`, `x.coint` or none.
  - Bartlett/Parzen bandwidth below 1 uses lag 1 with a negative weight.
- What the package does instead:
  - **FM-OLS** as cointReg, but the bias correction uses n = T − 1 (the
    observations in the sums) instead of T; ω̂_{u·v} and standard errors
    identical, coefficients differ by O(1/T) (< 2e-3 in the tests).
  - **DOLS**: Δx_{t−j}, j = −leads..lags, Δx_t always included; selection
    by AIC/BIC over 0..kmax (kmax = ⌊4(T/100)^{1/4}⌋) with every pair on the
    common sample t = 2 + kmax..T − kmax and all coefficients counted; the
    chosen pair is refitted on its maximal sample. Covariance: long-run
    variance of the residuals × (Z′Z)⁻¹ (as cointReg). Identical to cointReg
    for given leads/lags ≥ 1 and to the regression written out with `lm`.
  - **IM-OLS** (type 1): identical to cointReg, with the covariance
    ω̂_{u·v}(Z′Z)⁻¹C′C(Z′Z)⁻¹ correctly scaled. ω̂_{u·v} from the FM-OLS
    long-run quantities of the OLS residuals, as cointReg: kernel
    inference, **not fixed-b**. Fixed-b critical values: later, by our own
    simulation **[open]**. **[check]** whether Vogelsang and Wagner (2014)
    recommend ω̂_{u·v} from the IM-OLS residuals instead.
  - Kernels Bartlett (default), Parzen, QS; bandwidth Andrews (AR(1) per
    column, default), Newey–West or a number. p-values standard normal.
  - Coefficient names `(Intercept)`, `trend`, and the formula's terms.
  - `residuals()` is the equilibrium error ẑ_t = y_t − θ̂′(d_t, x_t) over the
    whole sample for every method (what the EG test and the ECM use).
- Monte Carlo (triangular DGP, corr(u, v) = 0.6, AR(0.5) error, T = 200,
  500 reps), coverage of the 95% interval for β: OLS 0.56, FM 0.88,
  DOLS 0.88, IM 0.91.
- Validation (tests): FM/D/IM against cointReg, DOLS against the regression
  written out with `lm`. Still to do: FM-OLS against EViews or MATLAB output
  for one dataset.

### 3.3 Engle–Granger test — `eg_test()`

```r
eg_test(y ~ x1 + x2, data, trend = c("c", "ct", "n"),
        lags = NULL, max_lags = NULL, criterion = "AIC",
        normalisations = FALSE)
```

- Step 1 = `coint_reg(method = "OLS")`; Step 2 = ADF with **no deterministic
  terms** on ẑ_t, through `unitroottools::adf(trend = "n")` so lag selection,
  fixed-sample IC comparison and residual checks are the ones students already
  know from ADF_T. Only the reference distribution changes.
- **Critical values and p-values depend on n and on Step 1's trend**:
  - Critical values: **MacKinnon (2010), Table 2**, response surfaces
    c(N) = β∞ + β₁/T + β₂/T² (+ β₃/T³), N = 1…12, cases "nc" (N = 1 only),
    "c", "ct", "ctt". Stored in `sysdata.rda` exactly as
    `mackinnon_tau` in unitroottools; N = 1 must reproduce unitroottools'
    values (unit test).
  - p-values: **own simulated tables, no urca** (decided 25 Nov 2026).
    Quantiles of the residual τ statistic on a grid of T, n = 1…12 and
    Step 1 trend, simulated with the Rcpp set-up of unitroottools'
    `data-raw/` (as for rho, phi, GLS); response surfaces per quantile;
    p-value by MacKinnon's (1996) local regression of Φ⁻¹(p) on the
    statistic. n = 1 is the plain ADF τ, so the same table can replace
    `urca::punitroot()` in unitroottools, which is urca's only real use
    there (`tau_p_value()`; other mentions are comments). Otherwise urca
    stays installed through unitroottools. Validation: critical values
    against MacKinnon (2010) Table 2 and p-values against `urca::punitroot()`
    (n = 1) and statsmodels' `coint()` (Suggests/tests only).
    The same machinery produces the ECM-test tables (§3.5).
    (Not used: urca's internal `.urcval()` with the MacKinnon 1996 tables
    for niv = 1…12.)
  - Printed with a line saying which table was used and that the plain ADF
    critical values do **not** apply (the chapter's main warning).
- `normalisations = TRUE`: runs the test with each variable as the dependent
  one and prints them side by side (EG's dependence on normalisation, the
  chapter's argument for Johansen). Cheap and very useful in class.
- **Phillips–Ouliaris** Z_t, Z_α later (§9): same Step 1, Phillips–Perron
  statistics from unitroottools on ẑ_t; Z_t uses the same MacKinnon tables,
  Z_α the `itt = 2` tables.

### 3.4 Cointegrated ECM (Engle–Granger two-step) — `coint_ecm()`

```r
coint_ecm(fit, p = 1, q = 0, data = NULL, controls = NULL)
# fit: a coint_reg object (any method); ẑ_{t-1} from it
```

Δy_t = δ₀ + b ẑ_{t−1} + Σ γ_j Δy_{t−j} + Σ δ_ij Δx_{i,t−j} + controls + ε_t.

- **Returns an adltools `adl` object** (Δy as dependent, Δx with their lags,
  ẑ_{t−1} as an unlagged regressor after the bar, `y ~ x | X` syntax). Then
  every ADL3/ADL4 tool works unchanged, which is exactly what the chapter
  says: `adl_diagnostics()`, `cusum_test()`, `cusumsq_test()`,
  `chow_break_test()`, `chow_forecast_test()`, `hansen_test()`,
  `adl_outliers()`, `adl_gets()` for the short-run lags. No duplicated code.
  Needs a check that adltools' "structural" functions refuse this object
  cleanly (adltools §3.5 says they should).
- δ₀ included by default (chapter's practitioner note); `intercept = FALSE`
  possible.
- Printed with b, its standard error and the statement that its t-ratio here
  is conventional **because ẑ is I(0) under cointegration**; it is not the
  ECM cointegration test (§3.5).
- Weak exogeneity of x for β: Δx_i on ẑ_{t−1} and lagged differences, test
  κ = 0. adltools' `weak_exogeneity_test()` builds ec from k̂; here ẑ comes
  from `coint_reg`. **[open]**: an `ec =` argument added to adltools'
  function, or a `coint_weak_exogeneity_test()` here. Proposal: the
  argument in adltools (one engine).
- Dynamic forecasting from the cointegrated ECM (chapter, "Forecasting with
  the ECM"): v0.2+ through a levels-reconstruction wrapper over
  `adl_forecast()`, giving ŷ in levels and the comparison with a
  differences-only model. Can wait.

### 3.5 ECM-based cointegration test — `ecm_test()`

```r
ecm_test(y ~ x1 + x2, data, trend = c("c", "ct"), p = 1, q = 0,
         form = c("unrestricted", "unit"))
```

- One-step (single-equation) test of Kremers, Ericsson and Dolado (1992),
  Banerjee, Dolado and Mestre (1998), Ericsson and MacKinnon (2002):
  - `form = "unrestricted"`: Δy_t = δ₀ + θ y_{t−1} + Σπ_i x_{i,t−1} +
    short-run terms; t-ratio on θ.
  - `form = "unit"`: the chapter's "transparent" form, θ(y − x)_{t−1} + π x_{t−1}
    (n = 2; identical t-ratio on θ to the unrestricted form for the same
    regressors — unit test).
- **Plain OLS standard errors, never HAC** (chapter's callout): no `cov_type`
  argument at all, so it cannot be done wrong.
- Distribution: Ericsson and MacKinnon (2002) response surfaces in (n, T,
  trend, number of weakly exogenous regressors k). **[open]**: transcribe
  their coefficients (must be checked line by line against the paper; their
  `ecmtest` spreadsheet is the other published source) or simulate them
  with the unitroottools Rcpp set-up and check against the paper's
  tables. The chapter's Table 2 should then be regenerated from whichever
  source is adopted.
- Must not be confused with the two-step t-ratio on ẑ_{t−1} in `coint_ecm()`;
  each print method says which one it is.

### 3.6 Teaching simulators — `simulate_coint()`, `simulate_eg()`

- `simulate_coint(T, n, beta, rho_z, endogeneity, ...)`: DGPs of the chapter
  in one place — independent random walks, common-trend system, AR(1)
  equilibrium error, triangular Phillips system with correlated innovations
  (to make OLS t-ratios fail and FM/DOLS work), slope break (for GH).
- `simulate_eg(T, M, n, tests = c("adf", "eg", "ecm"))`: the chapter's
  three-test Monte Carlo, with `print`/`plot` like unitroottools'
  `simulate_df()`. Replaces the long code chunks in the chapter with one call
  while the code stays visible in the package source.
- Spurious-regression simulator (`simulate_spurious(T, M)`): R², |t| > 1.96
  share, growth of the rejection rate with T.

## 4. Test results and printing

Same class layout as adltools/forecastingtools (`statistic`, `p.value`,
`distribution`, `critical_values`, `null_hypothesis`, `reject_null`,
`alpha`, `nobs`, plus `n_vars`, `trend`, `lags`, `table_source`).
Class name `cointtools_test` **[open]**, or reuse `adltools_test` by
importing its constructor (would need exporting it from adltools).
`tidy()`/`glance()` via generics.

## 5. Conventions shared with the suite

- snake_case; tests end in `_test`; prefix `coint_` where a name would clash.
- Input types and NA handling as adltools §12; `start`/`end` as in dynlm.
- Tables are data frames; plots ggplot2/patchwork (`interactive = TRUE` for
  plotly) or base, matching what adltools ended up using **[check]**.
- MIT license; DESCRIPTION author block with ORCID and affiliation as in
  adltools.
- `Remotes: maarranz/UnitRoot-Tools-R, maarranz/ADL-Tools-R` while the
  siblings are GitHub-only.

## 6. Dependencies

Imports now: stats, tools, zoo. To add with the functions that need them:
unitroottools, adltools, generics (no urca, no cointReg; §3.2–3.3).
Suggests: cointReg (reference in tests), testthat, knitr, rmarkdown, dynlm (cross-checks),
strucchange, cointmonitoR (§9), tibble, tsibble, xts, fredr (vignette
only, not run on check).

## 7. Data

- `term_structure`: monthly TB3MS and GS10, snapshot from FRED (vignettes and
  tests cannot call FRED). Same series as the chapter's applied example.
- `price_indices`: log PCEPI and CPI snapshot for the CI(2,1) companion
  file — only if §10 is taken up.
- A simulated dataset with known α, β and b (like `adl_sim`).
- Possibly the Hendry–Ericsson UK narrow-money data, linking the
  money-demand replication **[open]**.

## 8. Validation plan

- MacKinnon (2010) critical values: reproduce Table 2 to all decimals; N = 1
  equals unitroottools.
- EG statistic: equals `ur.df(residuals(lm(...)), type = "none")` for the
  same lags; p-values against statsmodels' `coint()` (MacKinnon 1994
  asymptotic; documented differences) and EViews.
- ECM test: size by simulation (n = 2, 3; T = 100, 500) against the adopted
  tables.
- coint_reg: done for FM/D/IM against cointReg and DOLS against `lm`;
  coverage in the triangular DGP recorded in §3.2.
- Development in the Linux workspace; `R CMD check --as-cran` on the M1 Pro
  and CI, as for adltools.

## 9. Roadmap (proposal)

- **v0.1.0** — input handling (from adltools), `integration_order()`,
  `coint_reg()` (OLS/FM/D/IM, own code) and `coint_wald_test()` — DONE,
  `eg_test()` with MacKinnon (2010) critical values and simulated
  p-values (§3.3), `normalisations`, `simulate_coint()`, `simulate_spurious()`,
  `term_structure` data, a vignette following the chapter's applied example
  up to Step 3.
- **v0.2.0** — `coint_ecm()` returning an `adl` object (diagnostics and
  stability from adltools), weak exogeneity via adltools, `ecm_test()`,
  `simulate_eg()`, vignette extended to Steps 4–5.
- **v0.3.0** — Gregory–Hansen (C, C/T, C/S; ADF*, Z_t*, Z_α*; their 1996
  Table 1 critical values; grid plot of ADF(τ)); Phillips–Ouliaris.
- **Later** — ECM levels forecasts; Hatemi-J two-break test; monitoring
  through `cointmonitoR` (Suggests; the chapter's "Going further");
  bounds test **[open]**.
- **Not planned** — Johansen/VECM (`urca::ca.jo()`).

## 9b. Order of work across the suite

1. **unitroottools — DONE (25 Nov 2026, decision D21).** Own τ1/τ2/τ3
   surfaces by T and k; urca only in Suggests. `data-raw/tau_kernel.cpp`
   and `tau_simulate.R` take `nvars`: with nvars > 1 they give the
   Engle–Granger residual statistics (Step 1 over t = 1..T with cases
   n/c/ct, residual ADF without deterministics over t = k+2..T, maximal
   sample as `adf()`). cointegrationtools reuses them: one run per n,
   output with an `nvars` column, surfaces fitted as for τ. No ctt case.
2. **adltools**: `ec =` argument in `weak_exogeneity_test()`.
3. **cointegrationtools v0.1.0** (needs 1).
4. **cointegrationtools v0.2.0**: ECM-test tables (same machinery as 1),
   `coint_ecm()`, `ecm_test()`.
5. Bounds test in adltools, re-exported here.

Simulation design to fix before step 1: grid of T, number of replications,
trend cases (n, c, ct; ctt?), and the p-value method near the tails.

## 10. Points to settle in the chapter as the package is written

- The chapter calls `ur.df(..., selectlags = "AIC")` on residuals and then
  compares with EG tables by hand; `eg_test()` replaces that and gives the
  p-value the chapter says R lacks.
- The chapter states no software returns ECM-test p-values; `ecm_test()`
  would change that sentence.
- `get_theta_x()`/`get_theta_const()` disappear once `coint_reg()` exists.
- CI(2,1) (PCEPI/CPI): needs I(1)-level tables applied to differences or
  ad hoc simulation; keep out of the package unless the companion file
  needs a function.

## 11. Open questions

1. ~~EG p-values~~: decided, own simulated tables, no urca (§3.3). The
   tables are split: unitroottools ships only n = 1 (plain ADF τ, replacing
   `urca::punitroot()`); the n ≥ 2 Engle–Granger tables live only in
   cointegrationtools. The Rcpp kernel is data-raw code (not shipped), so
   the same script can sit in both repos; cointegrationtools also runs
   n = 1 as a check against unitroottools.
2. ~~ECM-test distribution~~: decided (25 Nov), own tables. Simulated once
   in `data-raw/`, response-surface coefficients stored in `sysdata.rda`;
   a test only evaluates the surface. Checked against Ericsson–MacKinnon
   (2002).
3. ~~Weak exogeneity~~: decided, available in both packages. One engine:
   `ec =` argument in adltools' `weak_exogeneity_test()`; the
   cointegrationtools function passes ẑ from `coint_reg()` to it.
4. Test class: undecided. Leaning to an own `cointtools_test` with the same
   fields, so no adltools internals are needed.
5. Bounds test: in both packages, implemented once in adltools (the ARDL in
   levels) and re-exported here (`importFrom` + `export`). Own simulated
   tables, checked against Pesaran, Shin and Smith (2001) and Narayan
   (2005).
6. ~~Python~~: yes, later; not needed by students now. Export the simulated
   quantiles as `csv.gz` (as unitroottools already does) so both languages
   read one set of tables. Note: Python's `arch` (`arch.unitroot.cointegration`)
   has Engle–Granger, Phillips–Ouliaris, DOLS and FM-OLS, and statsmodels
   has `coint()`; no IM-OLS and no ECM test. Miguel uses both in his
   existing Python packages and dislikes them: the Python version should
   depend only on numpy/scipy (pandas for input), with its own estimators
   and the shared tables, not wrap `arch`/statsmodels.
7. No EViews dataset for now. Validation by simulated DGPs and cross-checks
   against cointReg, statsmodels/`arch` and published tables.

## 12. References (to be checked against publisher records)

Engle and Granger (1987); MacKinnon (1996, 2010); Kremers, Ericsson and
Dolado (1992); Banerjee, Dolado and Mestre (1998); Ericsson and MacKinnon
(2002); Phillips and Hansen (1990); Phillips and Loretan (1991); Saikkonen
(1991); Stock and Watson (1993); Vogelsang and Wagner (2014); Phillips and
Ouliaris (1990); Gregory and Hansen (1996); Wagner and Wied (cointmonitoR;
the chapter's key is `wagnerWied2015`, check the year against the journal
version); Aschersleben and Wagner (2016, cointReg). Keys as in
`referencesCointegration.bib`.
