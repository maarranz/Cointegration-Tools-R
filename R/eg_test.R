# Engle-Granger residual-based test -------------------------------------------------

#' Engle-Granger test for cointegration
#'
#' The two-step test of Engle and Granger (1987):
#'
#' 1. **Cointegrating regression.** OLS of the levels regression
#'    \eqn{y_t = \alpha + \beta' x_t + z_t} (with a trend if `trend = "ct"`,
#'    nothing if `trend = "n"`), as [coint_reg()] with `method = "OLS"`.
#' 2. **Unit-root test on the residuals.** ADF regression of
#'    \eqn{\Delta \hat z_t} on \eqn{\hat z_{t-1}} and \eqn{k} lagged
#'    differences, **without deterministic terms** (the Step 1 regression has
#'    already removed them), as `unitroottools::adf(trend = "n")`, with the
#'    the same lag-selection options (BIC by default).
#'
#' \eqn{H_0}: no cointegration (\eqn{\hat z_t} has a unit root);
#' \eqn{H_1}: cointegration. The statistic is the ADF t ratio, but its
#' distribution is not the ADF one: it depends on the number of variables
#' \eqn{n} in the cointegrating regression and on the deterministic terms of
#' Step 1, and is further to the left. Critical values and p-values come from
#' the package's own response surfaces (simulated by sample size, lag order,
#' \eqn{n} and Step 1 case; see `data-raw/`), not from the plain ADF tables.
#'
#' The sample size used by the surfaces is the length of a series whose
#' maximal ADF sample is the effective sample, \eqn{T = N + k + 1}, as in
#' unitroottools. Tables are available for \eqn{n = 2, \dots, 6}.
#'
#' The result depends on which variable is put on the left of Step 1. With
#' `normalisations = TRUE` the test is repeated with each variable as the
#' dependent one, and the results are shown side by side.
#'
#' @inheritParams coint_reg
#' @param lags Fixed number of lagged differences in the ADF regression;
#'   `NULL` selects it with `criterion`.
#' @param max_lags Largest order searched (default: Schwert's rule, as in
#'   unitroottools).
#' @param criterion Lag selection: `"BIC"` (default), `"AIC"`, `"MAIC"` or
#'   `"MBIC"`, as in `unitroottools::adf()`. The default differs from
#'   `adf()` (MAIC): on cointegrating residuals MAIC tends to choose long lags
#'   when the data are far from the null, which costs power (see
#'   `DESIGN.md`, section 3.3).
#' @param final_sample `"maximal"` (default) or `"common"`, as in
#'   `unitroottools::adf()`.
#' @param alpha Significance level for the decision.
#' @param normalisations If `TRUE`, also run the test with each variable as
#'   the dependent one.
#' @return A `"cointtools_test"` object with, besides the usual fields,
#'   `critical_values` (1%, 5%, 10%), `lags`, `nvars`, `trend`,
#'   `regression` (the Step 1 [coint_reg()] fit), `adf` (the Step 2
#'   unitroottools result: its own critical values and p-value are those of
#'   the plain ADF test and must not be used), `warnings` and, if requested,
#'   `normalisations` (a data frame).
#' @references
#' Engle, R. F. and Granger, C. W. J. (1987). Co-integration and error
#' correction: representation, estimation, and testing. *Econometrica* 55(2),
#' 251-276.
#'
#' MacKinnon, J. G. (1996). Numerical distribution functions for unit root
#' and cointegration tests. *Journal of Applied Econometrics* 11(6), 601-618.
#' @seealso [coint_reg()]
#' @examples
#' set.seed(1)
#' n <- 200
#' w <- cumsum(rnorm(n))
#' d <- data.frame(x = w + rnorm(n, sd = 0.4))
#' d$y <- 2 + 1.5 * d$x + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4))
#' eg_test(y ~ x, d)
#'
#' # two independent random walks: no cointegration
#' d$y2 <- cumsum(rnorm(n))
#' eg_test(y2 ~ x, d, normalisations = TRUE)
#' @export
eg_test <- function(formula, data, trend = c("c", "ct", "n"), lags = NULL,
                    max_lags = NULL, criterion = c("BIC", "AIC", "MAIC", "MBIC"),
                    final_sample = c("maximal", "common"), alpha = 0.05,
                    normalisations = FALSE, time = NULL, start = NULL, end = NULL) {
  trend <- match.arg(trend)
  criterion <- match.arg(criterion)
  final_sample <- match.arg(final_sample)
  res <- eg_core(formula, data, trend, lags, max_lags, criterion, final_sample, alpha,
                 time, start, end)
  if (isTRUE(normalisations)) {
    spec <- parse_coint_formula(formula)
    terms <- c(spec$lhs_label, spec$rhs_labels)
    res$normalisations <- do.call(rbind, lapply(seq_along(terms), function(i) {
      f <- stats::reformulate(terms[-i], response = str2lang(terms[i]),
                              env = environment(formula))
      r <- if (i == 1) res else suppressMessages(eg_core(f, data, trend, lags, max_lags,
                                        criterion, final_sample, alpha, time, start, end))
      data.frame(dependent = terms[i], statistic = r$statistic, lags = r$lags,
                 p.value = r$p.value, reject = r$reject_null)
    }))
  }
  res
}

#' One Engle-Granger test (no normalisation loop)
#' @noRd
eg_core <- function(formula, data, trend, lags, max_lags, criterion, final_sample,
                    alpha, time, start, end) {
  fit <- coint_reg(formula, data, method = "OLS", trend = trend, time = time,
                   start = start, end = end)
  nvars <- 1L + length(fit$xnames)
  if (!as.character(nvars) %in% names(eg_rs)) {
    stop("Engle-Granger tables are available for ",
         paste(range(as.integer(names(eg_rs))), collapse = " to "),
         " variables; this regression has ", nvars, ".", call. = FALSE)
  }
  z <- as.numeric(zoo::coredata(fit$residuals))
  a <- unitroottools::adf(z, trend = "n", lags = lags, max_lags = max_lags,
                          criterion = criterion, final_sample = final_sample,
                          level = alpha)
  k <- a$lags
  T_eq <- a$nobs + k + 1L
  rs <- eg_rs[[as.character(nvars)]][[trend]]
  cv <- rs_critical_values(rs, T_eq, k)
  p <- rs_p_value(rs, a$statistic, T_eq, k)
  # adf()'s own design warnings refer to the same (T, k) grid: unique() keeps one
  warn <- unique(c(a$warnings, rs_design_warnings(rs, T_eq, k), p_value_note(p)))
  sel <- a$lag_selection
  lag_txt <- if (is.null(sel)) paste0(k, " (given)") else
    paste0(k, " (", sel$method, ", max ", sel$kmax, ")")
  details <- c(
    `Step 1` = paste0("OLS of ", fit$yname, " on ", paste(fit$xnames, collapse = ", "),
                      "; deterministic: ", trend_labels[[trend]]),
    `Step 2` = paste0("ADF on the residuals, no deterministic terms; lags ", lag_txt),
    `Tables` = paste0("Engle-Granger, n = ", nvars, " variables, ", trend_labels[[trend]],
                      " in Step 1 (own response surfaces; not the ADF tables)"))
  new_coint_test(list(
    method = "Engle-Granger", null_hypothesis = "no cointegration (residuals have a unit root)",
    alternative = "cointegration (residuals stationary)",
    statistic = unname(a$statistic), distribution = "Engle-Granger tau",
    critical_values = cv, p.value = p, alpha = alpha,
    reject_null = isTRUE(p < alpha), nobs = a$nobs, lags = k, nvars = nvars,
    trend = trend, details = details, warnings = warn,
    regression = fit, adf = a, coefficients = NULL),
    title = "Engle-Granger cointegration test")
}
