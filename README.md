# cointegrationtools (R)

[![R-CMD-check](https://github.com/maarranz/Cointegration-Tools-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maarranz/Cointegration-Tools-R/actions/workflows/R-CMD-check.yaml)

Teaching-oriented R tools for **single-equation cointegration**, developed for **Applied Time Series Econometrics (ATSE)** at **Universidad Carlos III de Madrid**.

> **Status: version 0.1.0.** `integration_order()`, `coint_reg()` (OLS, FM-OLS, DOLS, IM-OLS), `coint_wald_test()`, `eg_test()`, the simulators and the term-structure vignette are available. The ECM functions (`coint_ecm()`, `ecm_test()`) are planned for 0.2.0; names and arguments may still change. The design decisions are recorded in [`DESIGN.md`](DESIGN.md).

## Scope

`cointegrationtools` covers cointegration as it is taught in the ATSE chapter on cointegration: one cointegrating relation among any number of I(1) variables, tested and estimated equation by equation. It reuses `unitroottools` for the unit-root steps and `adltools` for the diagnostics and stability of the error-correction model.

| Task | Planned functions |
|:--|:--|
| Integration order | `integration_order()`: unit-root tests on each series, in levels and differences (through `unitroottools`) |
| Cointegrating regression | `coint_reg()`: OLS, fully modified OLS, dynamic OLS and integrated modified OLS, one interface and one result class; `coint_wald_test()` for restrictions on the cointegrating vector |
| Engle–Granger test | `eg_test()`: residual-based ADF test with finite-sample critical values and p-values that depend on the number of variables and the deterministic terms; optional comparison across normalisations |
| ECM-based test | `ecm_test()`: the t-ratio on the error-correction term (Kremers, Ericsson and Dolado; Ericsson and MacKinnon), with its own critical values and p-values |
| Error-correction model | `coint_ecm()`: the two-step ECM, returned as an `adltools` model so that its diagnostics, CUSUM, Chow and Hansen tests apply directly; weak exogeneity of the regressors |
| Teaching simulations | `simulate_coint()`, `simulate_spurious()`, `simulate_eg()`: spurious regression, cointegrated systems, and the null distributions of the ADF, Engle–Granger and ECM statistics |

Critical values and p-values come from the package's own response surfaces, simulated once by sample size and lag order (code in `data-raw/`), as in `unitroottools`.

Later: Gregory–Hansen tests for a break in the cointegrating vector, Phillips–Ouliaris tests, and the bounds test (shared with `adltools`). Not planned: Johansen's system approach, available in `urca::ca.jo()`.

## Related packages

`cointegrationtools` is part of a suite of teaching packages for ATSE:

| Package | Repository | Covers |
|:--|:--|:--|
| arimatools | [ARIMA-Tools-R](https://github.com/maarranz/ARIMA-Tools-R) | ARIMA estimation, diagnostics, outliers and dynamic forecasting |
| forecastingtools | [Forecasting-Tools-R](https://github.com/maarranz/Forecasting-Tools-R) | Static forecasting, forecast accuracy, evaluation and comparison |
| unitroottools | [UnitRoot-Tools-R](https://github.com/maarranz/UnitRoot-Tools-R) | Unit-root and stationarity testing |
| adltools | [ADL-Tools-R](https://github.com/maarranz/ADL-Tools-R) | ADL models: estimation, inference, Granger causality, diagnostics, stability, exogeneity |
| cointegrationtools | Cointegration-Tools-R | Single-equation cointegration (this package) |

## Installation

```r
# install.packages("remotes")
remotes::install_github("maarranz/Cointegration-Tools-R", build_vignettes = TRUE)
library(cointegrationtools)
```

## Author

**Miguel A. Arranz**
Universidad Carlos III de Madrid
<maarranz@eco.uc3m.es>
ORCID: [0000-0002-5951-2284](https://orcid.org/0000-0002-5951-2284)

## Licence

MIT © Miguel A. Arranz
