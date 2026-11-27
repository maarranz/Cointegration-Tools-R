# Teaching simulations ----------------------------------------------------------------

#' Simulate cointegrated and non-cointegrated series
#'
#' The data-generating processes of the cointegration chapter, with every
#' series returned so that students can see what is behind each test:
#'
#' * `"cointegrated"`: \eqn{x_t} are independent random walks,
#'   \eqn{x_{it} = x_{i,t-1} + v_{it}}, and
#'   \eqn{y_t = \alpha + \beta' x_t + z_t} with a stationary AR(1)
#'   equilibrium error \eqn{z_t = \rho z_{t-1} + e_t}. `corr` makes
#'   \eqn{e_t} correlated with \eqn{v_{1t}}: the regressor is then
#'   endogenous, which is what makes OLS inference on \eqn{\beta} invalid and
#'   motivates FM-OLS, DOLS and IM-OLS.
#' * `"independent"`: \eqn{y_t} and \eqn{x_t} are independent random walks
#'   (no cointegration; the spurious-regression case).
#' * `"break"`: as `"cointegrated"`, but the slope changes from `beta` to
#'   `beta_post` after observation `floor(break_frac * T)` (cointegrated
#'   within each regime; the Gregory-Hansen case).
#'
#' All innovations are Gaussian with standard deviations `sd_x` (for
#' \eqn{v_t}) and `sd_z` (for \eqn{e_t}); random walks start at zero.
#'
#' @param T Number of observations.
#' @param beta Cointegrating coefficients, one per regressor (its length
#'   sets the number of regressors).
#' @param alpha Intercept of the long-run relation.
#' @param type `"cointegrated"`, `"independent"` or `"break"`.
#' @param rho AR(1) coefficient of the equilibrium error (\eqn{|\rho| < 1}).
#' @param corr Correlation between \eqn{e_t} and \eqn{v_{1t}}.
#' @param sd_x,sd_z Standard deviations of the innovations.
#' @param break_frac,beta_post Break fraction and post-break slope(s) for
#'   `type = "break"`.
#' @param seed Optional seed, for reproducible draws.
#' @return A data frame with `t`, `y`, `x` (or `x1`, `x2`, ...) and, except
#'   for `"independent"`, the true equilibrium error `z`.
#' @examples
#' d <- simulate_coint(200, beta = 1.5, rho = 0.5, corr = 0.6, seed = 1)
#' eg_test(y ~ x, d)
#' coint_reg(y ~ x, d, method = "FM")
#' @export
simulate_coint <- function(T = 200, beta = 1, alpha = 0,
                           type = c("cointegrated", "independent", "break"),
                           rho = 0.5, corr = 0, sd_x = 1, sd_z = 1,
                           break_frac = 0.5, beta_post = beta / 2, seed = NULL) {
  type <- match.arg(type)
  check_count(T, "T")
  if (T < 10) stop("`T` must be at least 10.", call. = FALSE)
  for (a in list(list(sd_x, "sd_x"), list(sd_z, "sd_z"))) {
    if (length(a[[1]]) != 1 || !is.numeric(a[[1]]) || !is.finite(a[[1]]) || a[[1]] <= 0) {
      stop("`", a[[2]], "` must be a single positive number.", call. = FALSE)
    }
  }
  if (!is.numeric(beta) || !length(beta) || any(!is.finite(beta))) {
    stop("`beta` must be a numeric vector of finite values.", call. = FALSE)
  }
  if (length(alpha) != 1 || !is.finite(alpha)) stop("`alpha` must be a single number.", call. = FALSE)
  if (length(rho) != 1 || !is.finite(rho) || length(corr) != 1 || !is.finite(corr)) {
    stop("`rho` and `corr` must be single numbers.", call. = FALSE)
  }
  if (abs(rho) >= 1) stop("`rho` must be between -1 and 1 (stationary error).", call. = FALSE)
  if (abs(corr) >= 1) stop("`corr` must be between -1 and 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  k <- length(beta)
  v <- matrix(stats::rnorm(T * k, sd = sd_x), T, k)
  e <- corr * v[, 1] / sd_x * sd_z + sqrt(1 - corr^2) * stats::rnorm(T, sd = sd_z)
  X <- apply(v, 2, cumsum)
  if (!is.matrix(X)) X <- matrix(X, T, k)
  colnames(X) <- if (k == 1) "x" else paste0("x", seq_len(k))
  out <- data.frame(t = seq_len(T), y = NA_real_, X)
  if (type == "independent") {
    out$y <- cumsum(stats::rnorm(T, sd = sd_x))
    return(out)
  }
  z <- as.numeric(stats::filter(e, rho, method = "recursive"))
  b <- matrix(beta, T, k, byrow = TRUE)
  if (type == "break") {
    if (length(beta_post) != k) stop("`beta_post` needs one value per regressor.", call. = FALSE)
    tb <- floor(break_frac * T)
    if (tb < 2 || tb > T - 2) stop("`break_frac` leaves no observations in a regime.", call. = FALSE)
    b[(tb + 1):T, ] <- matrix(beta_post, T - tb, k, byrow = TRUE)
  }
  out$y <- alpha + rowSums(b * X) + z
  out$z <- z
  out
}

#' Spurious regression by simulation
#'
#' Regresses one random walk on another, independent one, `M` times for each
#' sample size, and records \eqn{R^2} and the t ratio of the slope. With
#' I(1) variables and no cointegration the t ratio diverges: the share of
#' \eqn{|t| > 1.96} grows with \eqn{T} instead of staying at 5%
#' (Granger and Newbold, 1974; Phillips, 1986).
#'
#' @param T Sample sizes.
#' @param M Replications per sample size.
#' @param seed Optional seed.
#' @return A data frame of class `"simulate_spurious"` with, for each `T`,
#'   the mean and median \eqn{R^2}, the share of \eqn{|t| > 1.96}, and the
#'   median \eqn{|t|}; the draws are in `attr(x, "draws")`. `plot()` shows
#'   the rejection share against T.
#' @references
#' Granger, C. W. J. and Newbold, P. (1974). Spurious regressions in
#' econometrics. *Journal of Econometrics* 2(2), 111-120.
#'
#' Phillips, P. C. B. (1986). Understanding spurious regressions in
#' econometrics. *Journal of Econometrics* 33(3), 311-340.
#' @examples
#' s <- simulate_spurious(T = c(50, 100, 400), M = 500, seed = 1)
#' s
#' plot(s)
#' @export
simulate_spurious <- function(T = c(50, 100, 200, 500, 1000), M = 1000, seed = NULL) {
  check_count(M, "M")
  if (M < 1) stop("`M` must be at least 1.", call. = FALSE)
  if (!is.numeric(T) || !length(T) || any(!is.finite(T)) || any(T != round(T)) || any(T < 5)) {
    stop("`T` must be whole numbers of at least 5.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  draws <- do.call(rbind, lapply(T, function(n) {
    r <- t(vapply(seq_len(M), function(i) {
      y <- cumsum(stats::rnorm(n))
      x <- cumsum(stats::rnorm(n))
      xc <- x - mean(x)
      b <- sum(xc * y) / sum(xc^2)
      e <- y - mean(y) - b * xc
      s2 <- sum(e^2) / (n - 2)
      c(r2 = 1 - sum(e^2) / sum((y - mean(y))^2), t = b / sqrt(s2 / sum(xc^2)))
    }, numeric(2)))
    data.frame(T = n, r2 = r[, "r2"], t = r[, "t"])
  }))
  out <- do.call(rbind, lapply(split(draws, draws$T), function(d) {
    data.frame(T = d$T[1], mean_r2 = mean(d$r2), median_r2 = stats::median(d$r2),
               reject = mean(abs(d$t) > 1.96), median_abs_t = stats::median(abs(d$t)))
  }))
  rownames(out) <- NULL
  structure(out, class = c("simulate_spurious", "data.frame"), draws = draws, M = M)
}

#' @export
print.simulate_spurious <- function(x, digits = 3, ...) {
  title <- "Spurious regression: y on x, independent random walks"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat(attr(x, "M"), "replications per T. A valid test would reject 5% of the time.\n\n")
  tab <- as.data.frame(lapply(unclass(x), function(c) if (is.numeric(c)) signif(c, digits) else c))
  names(tab) <- c("T", "mean R2", "median R2", "share |t| > 1.96", "median |t|")
  print(tab, row.names = FALSE)
  invisible(x)
}

#' @export
plot.simulate_spurious <- function(x, ...) {
  graphics::plot(x$T, x$reject, type = "b", pch = 19, log = "x", ylim = c(0, 1),
                 xlab = "T (log scale)", ylab = "share of |t| > 1.96",
                 main = "Spurious regression: rejections of a true null", ...)
  graphics::abline(h = 0.05, lty = 2)
  invisible(x)
}

#' Null distributions of the ADF, Engle-Granger and ECM statistics
#'
#' The Monte Carlo of the cointegration chapter, in one call. Under the null
#' of no cointegration (`n` independent driftless random walks), computes in
#' each replication:
#'
#' * **ADF**: the Dickey-Fuller t ratio of the first regressor alone, with a
#'   constant;
#' * **EG**: the Engle-Granger statistic, the DF t ratio (no deterministic
#'   terms) of the residuals of y on a constant and the regressors;
#' * **ECM**: the t ratio of \eqn{\theta} in \eqn{\Delta y_t = \alpha + \theta
#'   y_{t-1} + \pi' x_{t-1} + \gamma' \Delta x_t + e_t}, with plain OLS
#'   standard errors (Kremers, Ericsson and Dolado, 1992; Ericsson and
#'   MacKinnon, 2002).
#'
#' All use `lags` lagged differences. The point is the comparison: all three
#' lie to the left of N(0, 1), EG more than ADF, and ECM in a place of its
#' own; none can be read against normal or ADF critical values.
#'
#' @param T Sample size.
#' @param M Replications.
#' @param n Number of variables (y and n - 1 regressors), 2 to 6.
#' @param lags Lagged differences in all three regressions.
#' @param seed Optional seed.
#' @return A data frame of class `"simulate_eg"` with the simulated 1%, 5% and
#'   10% quantiles of each statistic, the package's Engle-Granger critical
#'   values for the same case and N(0, 1) for reference; the draws are in
#'   `attr(x, "draws")`. `plot()` shows the three densities.
#' @examples
#' s <- simulate_eg(T = 200, M = 1000, seed = 1)
#' s
#' plot(s)
#' @export
simulate_eg <- function(T = 500, M = 2000, n = 2, lags = 0, seed = NULL) {
  check_count(T, "T")
  check_count(M, "M")
  check_count(lags, "lags")
  check_count(n, "n")
  if (n < 2 || n > 6) stop("`n` must be between 2 and 6.", call. = FALSE)
  if (T < 20 + 3 * lags + 2 * n) stop("`T` is too small for these lags.", call. = FALSE)
  if (M < 10) stop("`M` must be at least 10.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  tr <- function(Z, y, j) {
    f <- stats::lm.fit(Z, y)
    e <- f$residuals
    s2 <- sum(e^2) / (length(y) - ncol(Z))
    p <- ncol(Z)
    R <- f$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
    R[lower.tri(R)] <- 0
    V <- chol2inv(R)
    unname(f$coefficients[j] / sqrt(s2 * V[j, j]))
  }
  df_t <- function(v, const) {
    dv <- diff(v)
    N <- length(dv)
    rows <- (lags + 1):N
    Z <- cbind(if (const) 1, v[rows],
               if (lags) sapply(seq_len(lags), function(j) dv[rows - j]))
    tr(Z, dv[rows], 1 + const)
  }
  draws <- t(vapply(seq_len(M), function(i) {
    W <- apply(matrix(stats::rnorm(T * n), T, n), 2, cumsum)
    y <- W[, 1]; X <- W[, -1, drop = FALSE]
    adf <- df_t(X[, 1], TRUE)
    u <- stats::lm.fit(cbind(1, X), y)$residuals
    eg <- df_t(u, FALSE)
    dy <- diff(y); dX <- apply(X, 2, diff)
    if (!is.matrix(dX)) dX <- matrix(dX, ncol = n - 1)
    N <- length(dy)
    rows <- (lags + 1):N
    Z <- cbind(1, y[rows], X[rows, , drop = FALSE], dX[rows, , drop = FALSE],
               if (lags) sapply(seq_len(lags), function(j) dy[rows - j]),
               if (lags) do.call(cbind, lapply(seq_len(lags), function(j) dX[rows - j, , drop = FALSE])))
    ecm <- tr(Z, dy[rows], 2)
    c(ADF = adf, EG = eg, ECM = ecm)
  }, numeric(3)))
  pr <- c(0.01, 0.05, 0.10)
  rs <- eg_rs[[as.character(n)]][["c"]]
  tab <- rbind(
    t(apply(draws, 2, stats::quantile, probs = pr)),
    `EG tables` = unname(rs_critical_values(rs, T, lags)),
    `N(0,1)` = stats::qnorm(pr))
  colnames(tab) <- c("1%", "5%", "10%")
  out <- data.frame(statistic = rownames(tab), tab, check.names = FALSE, row.names = NULL)
  structure(out, class = c("simulate_eg", "data.frame"), draws = as.data.frame(draws),
            T = T, M = M, n = n, lags = lags)
}

#' @export
print.simulate_eg <- function(x, digits = 3, ...) {
  title <- "Null distributions: ADF, Engle-Granger and ECM statistics"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("T = ", attr(x, "T"), ", ", attr(x, "M"), " replications, n = ", attr(x, "n"),
      " independent random walks, ", attr(x, "lags"), " lagged difference(s).\n\n", sep = "")
  tab <- as.data.frame(unclass(x), check.names = FALSE)
  tab[-1] <- lapply(tab[-1], round, digits)
  print(tab, row.names = FALSE)
  cat("\nRows 1-3: simulated quantiles. 'EG tables': the package's critical values",
      "\nfor the same case (compare with row EG). None of the three is N(0, 1).\n")
  invisible(x)
}

#' @export
plot.simulate_eg <- function(x, ...) {
  d <- attr(x, "draws")
  dens <- lapply(d, stats::density)
  xr <- range(unlist(lapply(dens, `[[`, "x")), -4, 4)
  yr <- c(0, max(unlist(lapply(dens, `[[`, "y")), stats::dnorm(0)))
  cols <- c(ADF = "#0f6e56", EG = "#378add", ECM = "#993c1d")
  graphics::plot(NA, xlim = xr, ylim = yr, xlab = "statistic", ylab = "density",
                 main = "ADF vs Engle-Granger vs ECM under the null", ...)
  for (s in names(dens)) graphics::lines(dens[[s]], col = cols[[s]], lwd = 2)
  graphics::curve(stats::dnorm(x), add = TRUE, lty = 2, col = "grey40", lwd = 2)
  graphics::legend("topleft", c("ADF (single series)", "EG (residuals)", "ECM t ratio", "N(0,1)"),
                   col = c(cols, "grey40"), lwd = 2, lty = c(1, 1, 1, 2), bty = "n")
  invisible(x)
}
