# cointegrationtools 0.1.0

* `integration_order()`: order of integration of each variable from
  `unitroottools::ur_strategy()` on the differences and the levels, read from
  the highest order down (Dickey and Pantula, 1987); BIC lag selection.
* `simulate_coint()`: the chapter's data-generating processes (cointegrated
  with AR(1) error and optional endogeneity, independent random walks, slope
  break). `simulate_spurious()`: the spurious-regression Monte Carlo, with a
  `plot()` method.
* `term_structure`: monthly TB3MS and GS10 from FRED (`data-raw/term_structure.R`).
* Vignette "Cointegration in the term structure of interest rates".

* `coint_reg()`: cointegrating regression by static OLS, fully modified OLS
  (Phillips and Hansen, 1990), dynamic OLS (Saikkonen, 1991; Stock and
  Watson, 1993) and integrated modified OLS (Vogelsang and Wagner, 2014),
  with one formula interface and one result class (`print()`, `summary()`,
  `coef()`, `vcov()`, `residuals()`, `fitted()`, `nobs()`). Own
  implementation, checked against cointReg. DOLS chooses leads and lags on a
  common sample. Inputs: ts, zoo/xts, data frames, tibbles, tsibbles.
* `eg_test()`: Engle-Granger two-step test. Step 2 is `unitroottools::adf()`
  without deterministic terms on the Step 1 residuals, lags chosen by BIC by
  default (MAIC, `adf()`'s default, costs power on cointegrating residuals);
  critical values and p-values from the package's own response surfaces for
  2 to 6 variables and Step 1 with none, a constant or a constant and trend
  (simulated in `data-raw/`, checked against MacKinnon 1996 and by size in
  fresh simulations). `normalisations = TRUE` repeats the test with each
  variable as the dependent one.
* `coint_wald_test()`: Wald tests of linear restrictions on the cointegrating
  vector, written as text (`"x = 1"`) or as `R` and `r`.
* Package skeleton: DESCRIPTION, help page, tests, R CMD check on GitHub
  Actions.
