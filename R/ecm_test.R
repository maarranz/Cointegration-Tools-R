# ECM-based cointegration test -----------------------------------------------------

#' Error-correction test for cointegration
#'
#' The single-equation test of Kremers, Ericsson and Dolado (1992), Banerjee,
#' Dolado and Mestre (1998) and Ericsson and MacKinnon (2002). Estimate the
#' unrestricted error-correction model
#' \deqn{\Delta y_t = d_t + \theta y_{t-1} + \pi' x_{t-1} + \gamma_0' \Delta
#'   x_t + \sum_{j=1}^{p} (\phi_j \Delta y_{t-j} + \gamma_j' \Delta x_{t-j})
#'   + e_t}
#' by OLS and use the t ratio of \eqn{\theta}, with **plain OLS standard
#' errors** (not HAC: the tables are for that statistic).
#' \eqn{H_0: \theta = 0} (no error correction, no cointegration) against
#' \eqn{H_1: \theta < 0}. Writing the levels part as
#' \eqn{\theta (y_{t-1} - x_{t-1}) + (\pi + \theta) x_{t-1}}, as in the
#' chapter, gives the same t ratio for \eqn{\theta}.
#'
#' Critical values and p-values come from the package's own response
#' surfaces (simulated by sample size, `p`, number of variables and
#' deterministic case; see `data-raw/`), not from the ADF or Engle-Granger
#' tables. The test is usually more powerful than Engle-Granger, because it
#' does not impose the common-factor restriction implicit in the two-step
#' procedure.
#'
#' **Assumptions.** All variables are I(1) (see [integration_order()]), and
#' the regressors are weakly exogenous for the parameters of the equation: if
#' x also error-corrects, the single-equation statistic is not the right one
#' (Ericsson and MacKinnon, 2002). [coint_weak_exogeneity_test()] can be used
#' afterwards on a cointegrating fit.
#'
#' @inheritParams coint_reg
#' @param trend Deterministic terms: `"c"` constant (default) or `"ct"`
#'   constant and trend.
#' @param p Number of lagged differences of y and x. `NULL` (default) selects
#'   it by `criterion` over 0, ..., `max_p` on a common sample.
#' @param max_p Largest p searched (at most 4, the range of the tables).
#' @param criterion `"BIC"` (default) or `"AIC"`.
#' @param alpha Significance level for the decision.
#' @return A `"cointtools_test"` object with the usual fields,
#'   `critical_values`, `p`, `nvars`, `trend`, the ECM regression
#'   (`regression`, an `lm` object), the implied long-run coefficients
#'   \eqn{-\pi/\theta} (`long_run`), `assumption` and `warnings`.
#' @references
#' Banerjee, A., Dolado, J. J. and Mestre, R. (1998). Error-correction
#' mechanism tests for cointegration in a single-equation framework.
#' *Journal of Time Series Analysis* 19(3), 267-283.
#'
#' Ericsson, N. R. and MacKinnon, J. G. (2002). Distributions of error
#' correction tests for cointegration. *Econometrics Journal* 5(2), 285-318.
#'
#' Kremers, J. J. M., Ericsson, N. R. and Dolado, J. J. (1992). The power of
#' cointegration tests. *Oxford Bulletin of Economics and Statistics* 54(3),
#' 325-348.
#' @seealso [eg_test()], [coint_ecm()]
#' @examples
#' ecm_test(gs10 ~ tb3ms, term_structure)
#' d <- simulate_coint(200, type = "independent", seed = 1)
#' ecm_test(y ~ x, d)
#' @export
ecm_test <- function(formula, data, trend = c("c", "ct"), p = NULL, max_p = 4,
                     criterion = c("BIC", "AIC"), alpha = 0.05, time = NULL,
                     start = NULL, end = NULL) {
  trend <- match.arg(trend)
  criterion <- match.arg(criterion)
  check_prob(alpha, "alpha")
  check_count(p, "p")
  check_count(max_p, "max_p")
  ps <- if (is.null(p)) 0:max_p else p
  if (max(ps) > 4) {
    stop("The ECM-test tables cover p = 0 to 4 lagged differences.", call. = FALSE)
  }
  fit <- coint_reg(formula, data, method = "OLS", trend = trend, time = time,
                   start = start, end = end)
  nvars <- 1L + length(fit$xnames)
  if (!as.character(nvars) %in% names(ecm_rs)) {
    stop("ECM-test tables are available for ",
         paste(range(as.integer(names(ecm_rs))), collapse = " to "),
         " variables; this regression has ", nvars, ".", call. = FALSE)
  }
  y <- as.numeric(zoo::coredata(fit$y))
  X <- as.matrix(zoo::coredata(fit$x))
  T <- length(y)
  design <- function(p, rows) {
    dy <- c(NA, diff(y))
    dX <- rbind(NA, diff(X))
    Z <- cbind(`(Intercept)` = 1, trend = if (trend == "ct") rows,
               y_lag1 = y[rows - 1], X[rows - 1, , drop = FALSE], dX[rows, , drop = FALSE])
    k <- ncol(X)
    colnames(Z)[2 + (trend == "ct") + seq_len(k)] <- paste0(fit$xnames, "_lag1")
    colnames(Z)[2 + (trend == "ct") + k + seq_len(k)] <- paste0("d_", fit$xnames)
    for (j in seq_len(p)) {
      L <- cbind(dy[rows - j], dX[rows - j, , drop = FALSE])
      colnames(L) <- paste0("d_", c(fit$yname, fit$xnames), "_lag", j)
      Z <- cbind(Z, L)
    }
    list(Z = Z, dy = dy[rows])
  }
  ncoef <- function(p) 2 + (trend == "ct") + 2 * ncol(X) + p * nvars
  if (T - 1 - max(ps) - ncoef(max(ps)) < 10) {
    stop("Too few observations for the ECM test with p = ", max(ps), ".", call. = FALSE)
  }
  sel <- NULL
  if (length(ps) > 1) {
    rows <- seq.int(max(ps) + 2, T)
    n <- length(rows)
    pen <- if (criterion == "AIC") 2 else log(n)
    ic <- vapply(ps, function(q) {
      des <- design(q, rows)
      e <- stats::lm.fit(des$Z, des$dy)$residuals
      n * log(sum(e^2) / n) + pen * ncol(des$Z)
    }, numeric(1))
    p <- ps[which.min(ic)]
    sel <- data.frame(p = ps, value = ic)
  }
  rows <- seq.int(p + 2, T)
  des <- design(p, rows)
  Z <- des$Z
  reg <- stats::lm(des$dy ~ 0 + Z)
  names(reg$coefficients) <- colnames(Z)
  if (sqrt(mean(stats::residuals(reg)^2)) <= 1e-9 * max(abs(des$dy), .Machine$double.xmin)) {
    stop("The ECM regression fits exactly: there is nothing to test.", call. = FALSE)
  }
  b <- stats::coef(reg)
  se <- sqrt(diag(stats::vcov(reg)))
  stat <- b[["y_lag1"]] / se[["y_lag1"]]
  theta <- b[["y_lag1"]]
  pis <- b[paste0(fit$xnames, "_lag1")]
  long_run <- stats::setNames(-pis / theta, fit$xnames)
  nobs <- length(rows)
  T_eq <- nobs + p + 1L
  rs <- ecm_rs[[as.character(nvars)]][[trend]]
  cv <- rs_critical_values(rs, T_eq, p)
  pv <- rs_p_value(rs, stat, T_eq, p)
  bound <- min(rs$probs)
  reject <- if (pv <= bound + 1e-12 && alpha <= bound) NA else pv < alpha
  warn <- unique(c(rs_design_warnings(rs, T_eq, p), p_value_note(pv)))
  lag_txt <- if (is.null(sel)) paste0(p, " (given)") else
    paste0(p, " (", criterion, ", max ", max(ps), ")")
  details <- c(
    `Regression` = paste0("change in ", fit$yname, " on ", fit$yname, "(-1), ",
                          paste0(fit$xnames, "(-1)", collapse = ", "),
                          ", the changes in ", paste(fit$xnames, collapse = ", "),
                          "; deterministic: ", trend_labels[[trend]]),
    `Lags` = paste0(lag_txt, " lagged differences of all variables"),
    `Long run` = paste0(fit$yname, " = ", paste0(format(signif(long_run, 4)), " ", fit$xnames,
                                                 collapse = " + "), " (implied by -pi/theta)"),
    `Tables` = paste0("ECM test, n = ", nvars, " variables, ", trend_labels[[trend]],
                      " (own response surfaces; plain OLS t ratio)"))
  new_coint_test(list(
    method = "ECM", null_hypothesis = "no cointegration (theta = 0)",
    alternative = "cointegration (theta < 0)",
    statistic = unname(stat), distribution = "ECM t ratio", critical_values = cv,
    p.value = pv, alpha = alpha, reject_null = reject,
    p_bounded = pv <= bound + 1e-12 || pv >= 1 - bound - 1e-12,
    nobs = nobs, p = p, nvars = nvars, trend = trend, theta = theta,
    long_run = long_run, selection = sel, regression = reg, details = details,
    assumption = paste0("all ", nvars, " variables are I(1) and the regressors are ",
                        "weakly exogenous; check with integration_order()."),
    warnings = warn),
    title = "ECM test for cointegration")
}
