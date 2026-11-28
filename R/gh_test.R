# Gregory-Hansen test: cointegration with one break of unknown date ------------------

# Asymptotic critical values of ADF* and Z_t* (shared) and Z_alpha* at 1%, 5%
# and 10%, for m = 1..4 regressors: Gregory and Hansen (1996a), Table 1
# (models C, C/T, C/S), and Gregory and Hansen (1996b), Table 1 (C/S/T).
# Checked against COINT 0.0.4 (whose Z_alpha* values for m = 4 are
# overwritten by the C/S/T ones in models C, C/T and C/S, a bug) and against a
# simulation of the package (data-raw/gh_check.R).
gh_tables <- local({
  tab <- function(adf, za) {
    list(ADF = matrix(adf, 4, 3, byrow = TRUE, dimnames = list(1:4, c("0.01", "0.05", "0.1"))),
         Za = matrix(za, 4, 3, byrow = TRUE, dimnames = list(1:4, c("0.01", "0.05", "0.1"))))
  }
  list(
    `C` = tab(c(-5.13, -4.61, -4.34,  -5.44, -4.92, -4.69,  -5.77, -5.28, -5.02,  -6.05, -5.56, -5.31),
              c(-50.07, -40.48, -36.19,  -57.01, -46.98, -42.49,  -63.64, -53.58, -48.65,
                -70.18, -59.40, -54.38)),
    `C/T` = tab(c(-5.45, -4.99, -4.72,  -5.80, -5.29, -5.03,  -6.05, -5.57, -5.33,  -6.36, -5.83, -5.59),
                c(-57.28, -47.96, -43.22,  -64.77, -53.92, -48.94,  -70.27, -59.76, -54.94,
                  -76.95, -65.44, -60.12)),
    `C/S` = tab(c(-5.47, -4.95, -4.68,  -5.97, -5.50, -5.23,  -6.51, -6.00, -5.75,  -6.92, -6.41, -6.17),
                c(-57.17, -47.04, -41.85,  -68.21, -58.33, -52.85,  -80.15, -68.94, -63.42,
                  -90.35, -78.52, -72.56)),
    `C/S/T` = tab(c(-6.02, -5.50, -5.24,  -6.45, -5.96, -5.72,  -6.89, -6.32, -6.16,  -7.31, -6.84, -6.58),
                  c(-69.37, -58.58, -53.31,  -79.65, -68.43, -63.10,  -90.84, -78.87, -72.75,
                    -100.69, -88.47, -82.30)))
})

gh_model_labels <- c(`C` = "level shift", `C/T` = "level shift with trend",
                     `C/S` = "regime shift (intercept and slopes)",
                     `C/S/T` = "regime and trend shift")

#' Gregory-Hansen test for cointegration with a structural break
#'
#' The residual-based tests of Gregory and Hansen (1996a, b) of the null of no
#' cointegration against cointegration with one break, at an unknown date, in
#' the cointegrating regression. With \eqn{D_t = 1(t > t_b)}:
#'
#' * `"C"`, level shift: \eqn{y_t = \mu_1 + \mu_2 D_t + \beta' x_t + z_t};
#' * `"C/T"`, level shift with trend: adds \eqn{\delta t};
#' * `"C/S"`, regime shift: \eqn{y_t = \mu_1 + \mu_2 D_t + \beta_1' x_t +
#'   \beta_2' x_t D_t + z_t};
#' * `"C/S/T"`, regime and trend shift: adds \eqn{\delta_1 t + \delta_2 t D_t}
#'   to `"C/S"` (Gregory and Hansen, 1996b).
#'
#' For every break date \eqn{t_b} in the middle of the sample (a fraction
#' `trim` is left out at each end) the model is fitted by OLS and three
#' statistics are computed on its residuals: the ADF t ratio, as in
#' [eg_test()] (no deterministic terms, lags selected at each date), and the
#' Phillips-Ouliaris \eqn{Z_t} and \eqn{Z_\alpha}, as in [po_test()]. The
#' test statistics are the smallest values over the dates, \eqn{ADF^*},
#' \eqn{Z_t^*} and \eqn{Z_\alpha^*}; the date that gives \eqn{ADF^*} is the
#' estimated break.
#'
#' **Critical values are asymptotic**, from Gregory and Hansen's Table 1
#' (1%, 5% and 10%, up to four regressors); no p-value is computed, only the
#' interval in which it lies. Rejection is evidence of cointegration,
#' possibly with a break, not evidence of a break: a stable cointegrating
#' relation also rejects (and the Engle-Granger test, with smaller critical
#' values in absolute value, is then more powerful).
#'
#' @inheritParams eg_test
#' @inheritParams po_test
#' @param model `"C"`, `"C/T"`, `"C/S"` (default) or `"C/S/T"`.
#' @param statistic The statistic reported as the test: `"ADF"` (default),
#'   `"Zt"` or `"Za"`. All three are always computed.
#' @param trim Fraction of the sample excluded at each end from the break
#'   dates searched (default 0.15, as Gregory and Hansen).
#' @param max_lags Largest lag order searched at each date (default: Schwert's
#'   rule, as in unitroottools).
#' @param alpha Significance level for the decision: 0.01, 0.05 or 0.10, the
#'   levels of the tables.
#' @return A `"cointtools_test"` object (class `"coint_gh"` too) with, besides
#'   the usual fields, `statistics` (a data frame with the three statistics,
#'   their break dates and critical values), `break_date` (the first
#'   observation of the second regime at \eqn{ADF^*}, or at the reported
#'   statistic), `path` (the three statistics at every date searched) and
#'   `regression` (the OLS fit at the estimated break). `plot()` draws the
#'   path.
#' @references
#' Gregory, A. W. and Hansen, B. E. (1996a). Residual-based tests for
#' cointegration in models with regime shifts. *Journal of Econometrics*
#' 70(1), 99-126.
#'
#' Gregory, A. W. and Hansen, B. E. (1996b). Tests for cointegration in models
#' with regime and trend shifts. *Oxford Bulletin of Economics and Statistics*
#' 58(3), 555-560.
#' @seealso [eg_test()], [po_test()], [coint_recursive()]
#' @examples
#' set.seed(321)
#' n <- 200
#' w <- cumsum(rnorm(n))
#' d <- data.frame(x = w + rnorm(n, sd = 0.3))
#' d$y <- ifelse(seq_len(n) <= 100, 1, 0.4) * d$x + rnorm(n, sd = 0.3)
#' eg_test(y ~ x, d)
#' g <- gh_test(y ~ x, d)
#' g
#' plot(g)
#' @export
gh_test <- function(formula, data, model = c("C/S", "C", "C/T", "C/S/T"),
                    statistic = c("ADF", "Zt", "Za"), trim = 0.15, lags = NULL,
                    max_lags = NULL, criterion = c("BIC", "AIC", "MAIC", "MBIC"),
                    kernel = c("bartlett", "parzen", "qs"), bandwidth = "andrews",
                    alpha = 0.05, time = NULL, start = NULL, end = NULL) {
  model <- match.arg(model)
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  kernel <- match.arg(kernel)
  if (!isTRUE(alpha %in% c(0.01, 0.05, 0.10))) {
    stop("`alpha` must be 0.01, 0.05 or 0.10, the levels of the Gregory-Hansen tables.",
         call. = FALSE)
  }
  if (!is.numeric(trim) || length(trim) != 1 || !is.finite(trim) || trim <= 0 || trim >= 0.5) {
    stop("`trim` must be a number between 0 and 0.5.", call. = FALSE)
  }
  fit <- coint_reg(formula, data, method = "OLS", trend = "c", time = time,
                   start = start, end = end)
  m <- length(fit$xnames)
  if (m > 4) {
    stop("Gregory-Hansen tables are available for 1 to 4 regressors; this regression has ",
         m, ".", call. = FALSE)
  }
  y <- as.numeric(zoo::coredata(fit$y))
  X <- as.matrix(zoo::coredata(fit$x))
  colnames(X) <- fit$xnames
  T <- length(y)
  breaks <- seq.int(ceiling(trim * T), floor((1 - trim) * T))
  design <- function(tb) {
    D <- as.numeric(seq_len(T) > tb)
    tt <- seq_len(T)
    Z <- cbind(`(Intercept)` = 1, shift = D)
    if (model %in% c("C/T", "C/S/T")) Z <- cbind(Z, trend = tt)
    if (model == "C/S/T") Z <- cbind(Z, `trend:shift` = tt * D)
    Z <- cbind(Z, X)
    if (model %in% c("C/S", "C/S/T")) {
      XD <- X * D
      colnames(XD) <- paste0(fit$xnames, ":shift")
      Z <- cbind(Z, XD)
    }
    Z
  }
  if (T - ncol(design(breaks[1])) < 20 || length(breaks) < 2) {
    stop("Too few observations (", T, ") for the Gregory-Hansen test with model ",
         model, ".", call. = FALSE)
  }
  check_exact_fit(stats::lm.fit(design(breaks[1]), y)$residuals, y, fit$yname)
  # ADF at each date: a fast version of unitroottools::adf(trend = "n") for
  # BIC and AIC (same samples and criteria; checked in the tests), adf()
  # itself for MAIC and MBIC
  fast <- criterion %in% c("BIC", "AIC")
  if (fast && is.null(lags)) {
    kmax <- if (is.null(max_lags)) as.integer(floor(12 * (T / 100)^0.25)) else max_lags
    check_count(kmax, "max_lags")
  }
  path <- do.call(rbind, lapply(breaks, function(tb) {
    u <- stats::.lm.fit(design(tb), y)$residuals
    a <- if (fast) adf_n_fast(u, lags, kmax, criterion)
         else unitroottools::adf(u, trend = "n", lags = lags, max_lags = max_lags,
                                 criterion = criterion)
    ps <- po_stats(u, kernel, bandwidth)
    data.frame(tb = tb, ADF = unname(a$statistic), lags = a$lags, Zt = ps$Zt, Za = ps$Za)
  }))
  idx <- fit$index
  path$date <- idx[path$tb + 1L]
  tabs <- gh_tables[[model]]
  stats_tab <- do.call(rbind, lapply(c("ADF", "Zt", "Za"), function(s) {
    i <- which.min(path[[s]])
    cv <- (if (s == "Za") tabs$Za else tabs$ADF)[m, ]
    data.frame(statistic = s, value = path[[s]][i], break_date = as.character(path$date[i]),
               `cv 1%` = cv[[1]], `cv 5%` = cv[[2]], `cv 10%` = cv[[3]],
               p.value = gh_p_text(path[[s]][i], cv),
               reject = path[[s]][i] < cv[[format(alpha)]], check.names = FALSE)
  }))
  row <- stats_tab[stats_tab$statistic == statistic, ]
  ib <- which.min(path[[statistic]])
  tb <- path$tb[ib]
  Z <- design(tb)
  reg <- stats::lm(y ~ 0 + Z)
  names(reg$coefficients) <- colnames(Z)
  cv <- unlist(row[c("cv 1%", "cv 5%", "cv 10%")])
  names(cv) <- c("0.01", "0.05", "0.1")
  lag_txt <- if (is.null(lags)) paste0(criterion, " at each date, ", path$lags[ib],
                                         " at the break") else paste0(lags, " (given)")
  details <- c(
    `Model` = paste0(model, ", ", gh_model_labels[[model]], "; OLS of ", fit$yname, " on ",
                     paste(fit$xnames, collapse = ", ")),
    `Breaks` = paste0("dates from ", as.character(path$date[1]), " to ",
                      as.character(path$date[nrow(path)]), " (trim ", trim, "); ",
                      length(breaks), " regressions"),
    `Break` = paste0("second regime starts at ", as.character(path$date[ib]),
                     " (minimum of ", statistic, ")"),
    `ADF lags` = lag_txt,
    `Tables` = paste0("Gregory and Hansen (1996", if (model == "C/S/T") "b" else "a",
                      "), Table 1, m = ", m, if (m == 1) " regressor" else " regressors", ": asymptotic values"))
  out <- new_coint_test(list(
    method = "Gregory-Hansen",
    null_hypothesis = "no cointegration",
    alternative = paste0("cointegration with a break at an unknown date (", model, ")"),
    statistic = row$value, distribution = paste0("Gregory-Hansen ", statistic, "*"),
    critical_values = cv, p.value = NA_real_, p_text = row$p.value, alpha = alpha,
    reject_null = row$reject, nobs = T, nvars = m + 1L, model = model,
    statistic_name = statistic, statistics = stats_tab, break_date = path$date[ib],
    break_index = tb + 1L, path = path, trim = trim, regression = reg, details = details,
    assumption = paste0("all ", m + 1L, " variables are I(1); check with integration_order(). ",
                        "Rejection does not by itself show a break."),
    warnings = character()),
    title = "Gregory-Hansen test for cointegration with a break")
  class(out) <- c("coint_gh", class(out))
  out
}

#' ADF t ratio without deterministic terms, lags fixed or selected by BIC or
#' AIC on the common sample and the final regression on the maximal sample:
#' unitroottools::adf(trend = "n", final_sample = "maximal") without its
#' extras, for the many regressions of gh_test()
#' @noRd
adf_n_fast <- function(u, lags, kmax, criterion) {
  T <- length(u)
  du <- diff(u)
  build <- function(k, first) {
    t <- first:T
    X <- cbind(u[t - 1], vapply(seq_len(k), function(j) du[t - 1 - j], numeric(length(t))))
    list(X = X, y = du[t - 1])
  }
  if (is.null(lags)) {
    first <- kmax + 2L
    nobs <- T - first + 1L
    pen <- if (criterion == "BIC") log(nobs) else 2
    ic <- vapply(0:kmax, function(k) {
      d <- build(k, first)
      log(sum(stats::.lm.fit(d$X, d$y)$residuals^2) / nobs) + pen * k / nobs
    }, numeric(1))
    lags <- which.min(ic) - 1L
  }
  d <- build(lags, lags + 2L)
  f <- stats::.lm.fit(d$X, d$y)
  rss <- sum(f$residuals^2)
  v <- solve(crossprod(d$X))[1, 1]
  list(statistic = f$coefficients[[1]] / sqrt(rss / (nrow(d$X) - ncol(d$X)) * v), lags = lags)
}

#' The interval of an asymptotic p-value from three critical values
#' @noRd
gh_p_text <- function(stat, cv) {
  if (stat < cv[[1]]) "< 0.01" else if (stat < cv[[2]]) "0.01 to 0.05"
  else if (stat < cv[[3]]) "0.05 to 0.10" else "> 0.10"
}

#' Plot the Gregory-Hansen statistics over the break dates
#'
#' Draws the statistic at each candidate break date, with its 5% critical
#' value and the estimated break.
#'
#' @param x A result of [gh_test()].
#' @param statistic Which statistic to draw (default: the one reported).
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.coint_gh <- function(x, statistic = x$statistic_name, ...) {
  statistic <- match.arg(statistic, c("ADF", "Zt", "Za"))
  p <- x$path
  row <- x$statistics[x$statistics$statistic == statistic, ]
  i <- which.min(p[[statistic]])
  rng <- range(c(p[[statistic]], row$`cv 5%`))
  graphics::plot(p$date, p[[statistic]], type = "l", ylim = rng, xlab = "",
                 ylab = paste0(statistic, " at each break date"),
                 main = paste0("Gregory-Hansen, model ", x$model), ...)
  graphics::abline(h = row$`cv 5%`, lty = 2, col = "grey40")
  graphics::abline(v = p$date[i], lty = 3, col = "grey40")
  graphics::legend("topright", c("5% critical value", "estimated break"), lty = c(2, 3),
                   col = "grey40", bty = "n", cex = 0.8)
  invisible(x)
}
