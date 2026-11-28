# Phillips-Ouliaris residual-based tests ---------------------------------------------

#' Phillips-Ouliaris test for cointegration
#'
#' The residual-based tests of Phillips and Ouliaris (1990): the Phillips-Perron
#' statistics \eqn{Z_t} and \eqn{Z_\alpha} applied to the residuals of the
#' cointegrating regression. Instead of lagged differences (Engle-Granger),
#' the serial correlation of the residuals is handled by a kernel estimate of
#' the long-run variance. With \eqn{\hat u_t} the residuals of the
#' cointegrating regression, \eqn{\hat\rho} the OLS coefficient of
#' \eqn{\hat u_t} on \eqn{\hat u_{t-1}} (no deterministic terms) over
#' \eqn{N = T - 1} observations, \eqn{e_t} its residuals,
#' \eqn{s^2 = \sum e_t^2 / N} and \eqn{\lambda^2} the long-run variance of
#' \eqn{e_t},
#' \deqn{Z_\alpha = N(\hat\rho - 1) - \frac{\lambda^2 - s^2}{2 S}, \qquad
#'   Z_t = \frac{s}{\lambda} t_{\hat\rho} - \frac{\lambda^2 - s^2}{2 \lambda S^{1/2}},}
#' with \eqn{S = N^{-2} \sum \hat u_{t-1}^2} and \eqn{t_{\hat\rho}} the OLS
#' t ratio of \eqn{\hat\rho - 1}.
#'
#' \eqn{H_0}: no cointegration; both statistics reject for large negative
#' values. Asymptotically \eqn{Z_t} has the Engle-Granger distribution.
#' Critical values and p-values come from the package's own response surfaces
#' for each statistic, simulated (by \eqn{T}, \eqn{n} and the deterministic
#' case) with the default kernel and bandwidth; with other choices they are an
#' approximation.
#'
#' As for [eg_test()], **every variable must be I(1)**, and the result can
#' depend on which variable is put on the left.
#'
#' @inheritParams eg_test
#' @param statistic `"Zt"` (default) or `"Za"`.
#' @param kernel,bandwidth Long-run variance of the residual innovations, as
#'   in [coint_reg()] (default Bartlett with Andrews' bandwidth).
#' @return A `"cointtools_test"` object; `statistics` holds both statistics,
#'   and `rho`, `lambda2`, `s2` and `bandwidth` the pieces of the formulas.
#' @references
#' Phillips, P. C. B. and Ouliaris, S. (1990). Asymptotic properties of
#' residual based tests for cointegration. *Econometrica* 58(1), 165-193.
#' @seealso [eg_test()], [gh_test()]
#' @examples
#' po_test(gs10 ~ tb3ms, term_structure)
#' po_test(gs10 ~ tb3ms, term_structure, statistic = "Za")
#' @export
po_test <- function(formula, data, trend = c("c", "ct", "n"), statistic = c("Zt", "Za"),
                    kernel = c("bartlett", "parzen", "qs"), bandwidth = "andrews",
                    alpha = 0.05, time = NULL, start = NULL, end = NULL) {
  trend <- match.arg(trend)
  statistic <- match.arg(statistic)
  kernel <- match.arg(kernel)
  check_prob(alpha, "alpha")
  fit <- coint_reg(formula, data, method = "OLS", trend = trend, time = time,
                   start = start, end = end)
  nvars <- 1L + length(fit$xnames)
  if (!as.character(nvars) %in% names(po_rs)) {
    stop("Phillips-Ouliaris tables are available for ",
         paste(range(as.integer(names(po_rs))), collapse = " to "),
         " variables; this regression has ", nvars, ".", call. = FALSE)
  }
  u <- as.numeric(zoo::coredata(fit$residuals))
  y <- as.numeric(zoo::coredata(fit$y))
  check_exact_fit(u, y, fit$yname)
  ps <- po_stats(u, kernel, bandwidth)
  T <- length(u)
  rs <- po_rs[[as.character(nvars)]][[trend]][[statistic]]
  stat <- ps[[statistic]]
  cv <- rs_critical_values(rs, T, 0L)
  p <- rs_p_value(rs, stat, T, 0L)
  bound <- min(rs$probs)
  reject <- if (p <= bound + 1e-12 && alpha <= bound) NA else p < alpha
  default_lrv <- kernel == "bartlett" && identical(bandwidth, "andrews")
  warn <- unique(c(rs_design_warnings(rs, T, 0L), p_value_note(p),
                   if (!default_lrv) paste("the tables were simulated with the Bartlett kernel",
                                           "and Andrews' bandwidth; with other choices they are",
                                           "an approximation")))
  other <- setdiff(c("Zt", "Za"), statistic)
  details <- c(
    `Step 1` = paste0("OLS of ", fit$yname, " on ", paste(fit$xnames, collapse = ", "),
                      "; deterministic: ", trend_labels[[trend]]),
    `Step 2` = paste0("residuals on their lag, no deterministic terms; long-run variance: ",
                      kernel, " kernel, bandwidth ", format(signif(ps$bandwidth, 3)),
                      " (", if (is.numeric(bandwidth)) "fixed" else bandwidth, ")"),
    `Other` = paste0(other, " = ", format(signif(ps[[other]], 4))),
    `Tables` = paste0("Phillips-Ouliaris ", statistic, ", n = ", nvars, " variables, ",
                      trend_labels[[trend]], " in Step 1 (own response surfaces)"))
  new_coint_test(list(
    method = "Phillips-Ouliaris", null_hypothesis = "no cointegration (residuals have a unit root)",
    alternative = "cointegration (residuals stationary)",
    assumption = paste0("all ", nvars, " variables are I(1); check with integration_order()."),
    statistic = unname(stat), distribution = paste("Phillips-Ouliaris", statistic),
    critical_values = cv, p.value = p, alpha = alpha, reject_null = reject,
    p_bounded = p <= bound + 1e-12 || p >= 1 - bound - 1e-12,
    nobs = T - 1L, nvars = nvars, trend = trend, statistic_name = statistic,
    statistics = c(Zt = ps$Zt, Za = ps$Za), rho = ps$rho, lambda2 = ps$lambda2,
    s2 = ps$s2, bandwidth = ps$bandwidth, details = details, warnings = warn,
    regression = fit),
    title = "Phillips-Ouliaris cointegration test")
}

#' Phillips-Ouliaris Z_t and Z_alpha of a residual series (data-raw/po_kernel.cpp
#' computes the same with the default kernel and bandwidth; keep in step)
#' @noRd
po_stats <- function(u, kernel = "bartlett", bandwidth = "andrews") {
  T <- length(u)
  N <- T - 1L
  u1 <- u[-T]
  u0 <- u[-1]
  sxx <- sum(u1^2)
  rho <- sum(u0 * u1) / sxx
  e <- u0 - rho * u1
  s2 <- sum(e^2) / N
  sigma2 <- sum(e^2) / (N - 1)
  lr <- long_run_cov(matrix(e), kernel, bandwidth)
  lambda2 <- lr$Omega[1, 1]
  S <- sxx / N^2
  trho <- (rho - 1) / sqrt(sigma2 / sxx)
  list(Zt = sqrt(s2 / lambda2) * trho - (lambda2 - s2) / (2 * sqrt(lambda2 * S)),
       Za = N * (rho - 1) - (lambda2 - s2) / (2 * S),
       rho = rho, lambda2 = lambda2, s2 = s2, bandwidth = lr$bandwidth)
}

#' Stop when a cointegrating regression fits exactly
#' @noRd
check_exact_fit <- function(u, y, yname) {
  if (sqrt(mean(u^2)) <= 1e-9 * max(abs(y), .Machine$double.xmin)) {
    stop("The cointegrating regression fits exactly (residuals are zero up to ",
         "rounding): ", yname, " is an exact linear function of the regressors, ",
         "so there is no equilibrium error to test.", call. = FALSE)
  }
}
