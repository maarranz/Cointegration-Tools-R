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
