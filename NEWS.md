# cointegrationtools (development version)

* `coint_reg()`: cointegrating regression by static OLS, fully modified OLS
  (Phillips and Hansen, 1990), dynamic OLS (Saikkonen, 1991; Stock and
  Watson, 1993) and integrated modified OLS (Vogelsang and Wagner, 2014),
  with one formula interface and one result class (`print()`, `summary()`,
  `coef()`, `vcov()`, `residuals()`, `fitted()`, `nobs()`). Own
  implementation, checked against cointReg. DOLS chooses leads and lags on a
  common sample. Inputs: ts, zoo/xts, data frames, tibbles, tsibbles.
* `coint_wald_test()`: Wald tests of linear restrictions on the cointegrating
  vector, written as text (`"x = 1"`) or as `R` and `r`.
* Package skeleton: DESCRIPTION, help page, tests, R CMD check on GitHub
  Actions.
