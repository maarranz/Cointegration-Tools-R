# Error-correction model built on a cointegrating regression ----------------------

#' Error-correction model (Engle-Granger two-step)
#'
#' Second step of Engle and Granger (1987): with the equilibrium error
#' \eqn{\hat z_t} of a cointegrating regression ([coint_reg()], any method),
#' estimates by OLS
#' \deqn{\Delta y_t = \delta_0 + b\,\hat z_{t-1} + \sum_{i=1}^{p} \delta_i
#'   \Delta y_{t-i} + \sum_{j=0}^{q} \gamma_j' \Delta x_{t-j} + u_t .}
#' \eqn{b} is the speed of adjustment: the share of last period's
#' disequilibrium that is corrected this period (negative under
#' cointegration).
#'
#' The result is an `adltools` model (class `c("coint_ecm", "adl")`) for
#' \eqn{\Delta y_t} with \eqn{\hat z_{t-1}} as a regressor without lags
#' (named `ec`), so the diagnostics and stability tests of adltools apply
#' directly: `adltools::adl_diagnostics()`, `cusum_test()`, `cusumsq_test()`,
#' `chow_break_test()`, `hansen_test()`, `adl_outliers()`, and so on. The
#' long-run functions of adltools (`adl_long_run()`, `adl_dynamics()`)
#' describe this equation in differences and are not the cointegrating
#' relation, which is `fit`.
#'
#' Under cointegration \eqn{\hat z_t} is I(0), so the t ratio of \eqn{b}
#' has its usual distribution. This is **not** a test of cointegration:
#' without cointegration \eqn{\hat z_t} is I(1) and that t ratio is
#' nonstandard. Test cointegration first ([eg_test()]).
#'
#' @param fit A `"coint_reg"` object: the long-run relation.
#' @param p Number of lags of \eqn{\Delta y_t}.
#' @param q Largest lag of the \eqn{\Delta x_t} (from 0, i.e. \eqn{\Delta
#'   x_t} itself, to `q`).
#' @param intercept Include \eqn{\delta_0}? (Default `TRUE`: cheap insurance,
#'   see the chapter's practitioner note.)
#' @return An object of class `c("coint_ecm", "adl")`, with also `coint`
#'   (the long-run fit), `ec` (the equilibrium error \eqn{\hat z_t}, as a zoo
#'   series) and `adjustment` (estimate, standard error, t ratio and p-value
#'   of \eqn{b}).
#' @seealso [coint_reg()], [eg_test()], [coint_weak_exogeneity_test()]
#' @examples
#' lr <- coint_reg(gs10 ~ tb3ms, term_structure, method = "FM")
#' ecm <- coint_ecm(lr, p = 1, q = 1)
#' ecm
#' adltools::cusum_test(ecm)
#' @export
coint_ecm <- function(fit, p = 1, q = 0, intercept = TRUE) {
  if (!inherits(fit, "coint_reg")) stop("`fit` must be a coint_reg object.", call. = FALSE)
  check_count(p, "p")
  check_count(q, "q")
  dat <- ecm_data(fit)
  ylab <- dat$names[1]
  xlab <- dat$names[-1]
  lags <- c(stats::setNames(list(seq_len(p)), ylab),
            stats::setNames(rep(list(0:q), length(xlab)), xlab))
  f <- stats::as.formula(paste(ylab, "~", paste(xlab, collapse = " + "), "| ec"))
  out <- suppressMessages(adltools::adl(f, dat$z, lags = lags, intercept = intercept))
  b <- stats::coef(out)[["ec"]]
  se <- sqrt(stats::vcov(out)["ec", "ec"])
  out$coint <- fit
  out$ec <- fit$residuals
  out$adjustment <- c(estimate = b, std.error = se, statistic = b / se,
                      p.value = 2 * stats::pnorm(-abs(b / se)))
  out$ecm_order <- c(p = p, q = q)
  class(out) <- c("coint_ecm", class(out))
  out
}

#' Differences of y and x and the lagged equilibrium error, as a zoo matrix
#' @noRd
ecm_data <- function(fit) {
  labs <- c(fit$yname, fit$xnames)
  nm <- make.unique(paste0("d_", make.names(labs)))
  lv <- cbind(as.numeric(zoo::coredata(fit$y)), as.matrix(zoo::coredata(fit$x)))
  d <- rbind(NA, diff(lv))
  z <- as.numeric(zoo::coredata(fit$residuals))
  m <- cbind(d, c(NA, z[-length(z)]))
  colnames(m) <- c(nm, "ec")
  list(z = zoo::zoo(m[-1, , drop = FALSE], fit$index[-1]), names = nm)
}

#' @export
print.coint_ecm <- function(x, ...) {
  lr <- x$coint
  title <- "Error-correction model (Engle-Granger two-step)"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Long run:  ", lr$yname, " ~ ", paste(lr$xnames, collapse = " + "), " by ",
      method_labels[[lr$method]], " (", trend_labels[[lr$trend]], "); ec = its residual\n",
      sep = "")
  a <- x$adjustment
  cat("Adjustment: b = ", format(signif(a[["estimate"]], 4)), " (s.e. ",
      format(signif(a[["std.error"]], 3)), ", t = ", format(round(a[["statistic"]], 2)),
      "): about ", format(round(100 * abs(a[["estimate"]]), 1)),
      "% of a disequilibrium is corrected each period\n", sep = "")
  cat("The t ratio of b is standard only if there is cointegration; it is not a",
      "cointegration test.\n\n")
  NextMethod()
}

#' Weak exogeneity of a regressor for the cointegrating vector
#'
#' Tests whether a regressor \eqn{x} adjusts to the disequilibrium, with the
#' marginal regression (Urbain, 1992)
#' \deqn{\Delta x_t = c + \kappa\,\hat z_{t-1} + \sum_{j=1}^{L} (\phi_j
#'   \Delta x_{t-j} + \psi_j \Delta y_{t-j}) + e_t}
#' (with the lagged differences of the other regressors as well).
#' \eqn{H_0: \kappa = 0}: \eqn{x} does not error-correct, so it is weakly
#' exogenous for the cointegrating vector and single-equation analysis of
#' \eqn{y} loses nothing. Under cointegration \eqn{\hat z_{t-1}} is I(0) and
#' the z test is standard (asymptotically normal).
#'
#' @param fit A `"coint_reg"` or `"coint_ecm"` object.
#' @param x The regressor to test (default: the only one).
#' @param lags Number of lagged differences \eqn{L}.
#' @param alpha Significance level.
#' @return A `"cointtools_test"` object with `kappa` and the marginal
#'   regression (`marginal`, an `lm` object).
#' @references Urbain, J.-P. (1992). On weak exogeneity in error correction
#'   models. *Oxford Bulletin of Economics and Statistics* 54(2), 187-207.
#' @examples
#' lr <- coint_reg(gs10 ~ tb3ms, term_structure, method = "FM")
#' coint_weak_exogeneity_test(lr)
#' @export
coint_weak_exogeneity_test <- function(fit, x = NULL, lags = 1, alpha = 0.05) {
  if (inherits(fit, "coint_ecm")) fit <- fit$coint
  if (!inherits(fit, "coint_reg")) {
    stop("`fit` must be a coint_reg or coint_ecm object.", call. = FALSE)
  }
  check_count(lags, "lags")
  check_prob(alpha, "alpha")
  x <- x %||% if (length(fit$xnames) == 1) fit$xnames else
    stop("Several regressors: choose one with `x`.", call. = FALSE)
  if (!x %in% fit$xnames) stop("'", x, "' is not a regressor of the fit.", call. = FALSE)
  dat <- ecm_data(fit)
  m <- as.matrix(zoo::coredata(dat$z))
  N <- nrow(m)
  lagv <- function(v, j) c(rep(NA, j), v[seq_len(N - j)])
  dname <- dat$names[1 + match(x, fit$xnames)]
  cols <- list(dx = m[, dname], ec = m[, "ec"])
  for (v in dat$names) for (j in seq_len(lags)) {
    cols[[paste0(v, "_lag", j)]] <- lagv(m[, v], j)
  }
  df <- stats::na.omit(as.data.frame(cols))
  marg <- stats::lm(dx ~ ., data = df)
  b <- stats::coef(marg)[["ec"]]
  se <- sqrt(stats::vcov(marg)["ec", "ec"])
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  new_coint_test(list(
    method = "weak_exogeneity", null_hypothesis = paste0(x, " is weakly exogenous (kappa = 0)"),
    statistic = z, distribution = "normal", p.value = p, alpha = alpha,
    reject_null = p < alpha, nobs = nrow(df), kappa = b, marginal = marg,
    details = c(`Marginal` = paste0("change in ", x, " on the lagged equilibrium error and ",
                                    lags, " lag(s) of all differences"))),
    title = "Weak exogeneity test (regression for the change in x)")
}
