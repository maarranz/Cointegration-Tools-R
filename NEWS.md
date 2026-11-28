# cointegrationtools 0.2.0

* `ecm_test()`: the error-correction test for cointegration (Kremers,
  Ericsson and Dolado, 1992; Ericsson and MacKinnon, 2002), plain-OLS t ratio
  of y(-1) in the unrestricted ECM, with critical values and p-values from
  the package's own response surfaces for 2 to 6 variables, constant or
  constant and trend, and p = 0 to 4 lagged differences (BIC by default);
  reports the implied long-run coefficients.
* `simulate_eg()` also shows the ECM-test critical values.
* Vignette: Steps 4 (ECM test and ECM) and 5 (recursive estimates, CUSUM).

* `coint_ecm()`: the Engle-Granger two-step error-correction model on any
  `coint_reg()` fit, returned as an adltools model, so its diagnostics,
  CUSUM, Chow, Hansen and outlier tools apply directly; prints the speed of
  adjustment. `coint_weak_exogeneity_test()`: weak exogeneity of a regressor
  for the cointegrating vector (Urbain, 1992).
* `simulate_eg()`: the chapter's Monte Carlo of the ADF, Engle-Granger and
  ECM statistics under the null, with the package's EG critical values for
  comparison, and a density plot.
* `coint_recursive()`: recursive (growing-sample) or rolling estimates of the
  cointegrating regression with confidence bands, and a `plot()` method, for
  checking the stability of the long-run relation.

# cointegrationtools 0.1.1

Fixes from the audit of 26 November 2026 (details in `DESIGN.md`, section 13):

* `eg_test()` stops when the cointegrating regression fits exactly, states
  the I(1) assumption in its output, and gives no decision when the p-value
  is only a bound below `alpha`.
* `coint_reg(method = "D")`: when only `lags` or only `leads` is given, the
  other is selected with it held fixed; infeasible orders are refused.
* `coint_wald_test()`: text restrictions are parsed as linear expressions
  (never evaluated); named `R` columns are matched by name.
* Time columns given as text must be ISO dates and are checked for gaps.
* Formulas are evaluated in their own environment; regressors named
  `trend` or `(Intercept)` are refused; levels and simulator arguments are
  validated; p-values equal the level exactly at the tabulated critical
  values.

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
