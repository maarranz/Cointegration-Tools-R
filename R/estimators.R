# Estimators of the cointegrating regression ------------------------------------
#
#   y_t = delta' d_t + beta' x_t + u_t,   x_t = x_{t-1} + v_t,   t = 1, ..., T
#
# Internal functions on plain matrices: y (T vector), X (T x k, the I(1)
# regressors) and D (T x m deterministic terms, m may be 0). Each returns the
# coefficients theta = (delta', beta')', their covariance and method details.
# Written for this package and checked against cointReg (Aschersleben and
# Wagner, 2016); differences from cointReg are listed in DESIGN.md (section 3.2).

#' Least squares on a design matrix, with the inverse cross-product
#' @noRd
ls_fit <- function(Z, y) {
  qrZ <- qr(Z)
  if (qrZ$rank < ncol(Z)) {
    stop("The regressors of the cointegrating regression are collinear.", call. = FALSE)
  }
  coef <- qr.coef(qrZ, y)
  list(coef = coef, resid = as.numeric(y - Z %*% coef),
       ZtZinv = chol2inv(qr.R(qrZ)))
}

#' Static OLS
#' @noRd
coint_ols <- function(y, X, D) {
  Z <- cbind(D, X)
  f <- ls_fit(Z, y)
  s2 <- sum(f$resid^2) / (length(y) - ncol(Z))
  list(coef = f$coef, vcov = s2 * f$ZtZinv, nobs = length(y),
       details = list(sigma2 = s2))
}

#' Long-run quantities of (u_t, v_t) from the OLS residuals, t = 2, ..., T
#' @noRd
fm_long_run <- function(u_ols, X, kernel, bandwidth) {
  w <- cbind(u = u_ols[-1], diff(X))
  lr <- long_run_cov(w, kernel, bandwidth)
  Om <- lr_blocks(lr$Omega)
  De <- lr_blocks(lr$Delta)
  Ovv_inv_vu <- solve(Om$vv, Om$vu)
  list(lr = lr, Om = Om, De = De, Ovv_inv_vu = Ovv_inv_vu,
       omega_uv = as.numeric(Om$uu - Om$uv %*% Ovv_inv_vu))
}

#' Fully modified OLS (Phillips and Hansen, 1990)
#'
#' y+_t = y_t - dx_t' Omega_vv^{-1} Omega_vu and
#' theta = (Z'Z)^{-1} (Z'y+ - n [0; Delta+_vu]), Delta+_vu = Delta_vu -
#' Delta_vv Omega_vv^{-1} Omega_vu, over t = 2..T (n = T - 1 observations,
#' the number in the sums; cointReg uses T). Covariance omega_u.v (Z'Z)^{-1}.
#' @noRd
coint_fm <- function(y, X, D, kernel, bandwidth) {
  Z <- cbind(D, X)
  u_ols <- ls_fit(Z, y)$resid
  L <- fm_long_run(u_ols, X, kernel, bandwidth)
  dX <- diff(X)
  n <- nrow(dX)
  yplus <- y[-1] - as.numeric(dX %*% L$Ovv_inv_vu)
  Delta_plus <- L$De$vu - L$De$vv %*% L$Ovv_inv_vu
  Z1 <- Z[-1, , drop = FALSE]
  ZtZinv <- solve(crossprod(Z1))
  corr <- c(rep(0, ncol(D)), n * as.numeric(Delta_plus))
  coef <- as.numeric(ZtZinv %*% (crossprod(Z1, yplus) - corr))
  list(coef = coef, vcov = L$omega_uv * ZtZinv, nobs = n,
       details = list(omega_uv = L$omega_uv, kernel = kernel,
                      bandwidth = L$lr$bandwidth,
                      bandwidth_rule = L$lr$bandwidth_rule))
}

#' Design of the dynamic OLS regression for given lags and leads of dx
#'
#' Rows t whose dx_{t-j}, j = -leads..lags, all exist: t = 2 + lags, ...,
#' T - leads, restricted to `rows` when given (common sample).
#' @noRd
dols_design <- function(y, X, D, lags, leads, rows = NULL) {
  T <- length(y)
  dX <- rbind(NA, diff(X))                       # row t holds dx_t
  t_all <- rows %||% seq.int(2 + lags, T - leads)
  cols <- list(D[t_all, , drop = FALSE], X[t_all, , drop = FALSE])
  for (j in seq(-leads, lags)) {
    m <- dX[t_all - j, , drop = FALSE]
    colnames(m) <- paste0("d.", colnames(X),
                          if (j > 0) paste0(".lag", j) else if (j < 0) paste0(".lead", -j))
    cols[[length(cols) + 1]] <- m
  }
  list(Z = do.call(cbind, cols), y = y[t_all], rows = t_all)
}

#' Lags and leads of dx by information criterion on a common sample
#'
#' The free order(s) range over 0..kmax; an order given by the user is held
#' fixed. Every candidate is fitted on the same rows, t = 2 + max lags, ...,
#' T - max leads over the candidates, so that the criteria are comparable;
#' K counts all coefficients, deterministic terms included.
#' @noRd
dols_select <- function(y, X, D, kmax, ic, lags = NULL, leads = NULL) {
  T <- length(y)
  lag_grid <- if (is.null(lags)) 0:kmax else lags
  lead_grid <- if (is.null(leads)) 0:kmax else leads
  rows <- seq.int(2 + max(lag_grid), T - max(lead_grid))
  grid <- expand.grid(lags = lag_grid, leads = lead_grid)
  n <- length(rows)
  pen <- if (ic == "aic") 2 else log(n)
  grid$value <- apply(grid, 1, function(g) {
    des <- dols_design(y, X, D, g[["lags"]], g[["leads"]], rows)
    e <- ls_fit(des$Z, des$y)$resid
    n * log(sum(e^2) / n) + pen * ncol(des$Z)
  })
  best <- grid[which.min(grid$value), ]
  list(lags = best$lags, leads = best$leads, table = grid, nobs = n,
       ic = ic, kmax = kmax, fixed = c(lags = !is.null(lags), leads = !is.null(leads)))
}

#' Stop unless a DOLS design leaves enough residual degrees of freedom
#' @noRd
dols_feasible <- function(T, k, m, lags, leads, what) {
  n <- T - 1 - lags - leads
  p <- m + k * (lags + leads + 2)
  if (n - p < 10) {
    stop("DOLS with ", what, " needs more observations: ", n, " usable rows for ",
         p, " coefficients (at least 10 residual degrees of freedom are required). ",
         "Use fewer leads/lags or a smaller `kmax`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Dynamic OLS (Saikkonen, 1991; Stock and Watson, 1993)
#'
#' Regression of y_t on d_t, x_t and dx_{t-j}, j = -leads..lags (dx_t always
#' included), on its maximal sample. Covariance: long-run variance of the
#' residuals times (Z'Z)^{-1}, as in cointReg.
#' @noRd
coint_dols <- function(y, X, D, lags, leads, kmax, ic, kernel, bandwidth) {
  sel <- NULL
  T <- length(y)
  if (is.null(lags) || is.null(leads)) {
    kmax <- kmax %||% floor(4 * (T / 100)^(1 / 4))
    L <- lags %||% kmax
    Q <- leads %||% kmax
    dols_feasible(T, ncol(X), ncol(D), L, Q,
                  paste0("lags up to ", L, " and leads up to ", Q, " (selection)"))
    sel <- dols_select(y, X, D, kmax, ic, lags, leads)
    lags <- sel$lags
    leads <- sel$leads
  }
  dols_feasible(T, ncol(X), ncol(D), lags, leads, paste0(lags, " lags and ", leads, " leads"))
  des <- dols_design(y, X, D, lags, leads)
  f <- ls_fit(des$Z, des$y)
  lr <- long_run_cov(matrix(f$resid, ncol = 1, dimnames = list(NULL, "u")),
                     kernel, bandwidth)
  omega <- as.numeric(lr$Omega)
  p <- ncol(D) + ncol(X)
  V <- omega * f$ZtZinv
  list(coef = f$coef[1:p], vcov = V[1:p, 1:p, drop = FALSE], nobs = length(des$y),
       details = list(omega_u = omega, kernel = kernel, bandwidth = lr$bandwidth,
                      bandwidth_rule = lr$bandwidth_rule, lags = lags,
                      leads = leads, selection = sel,
                      dx_coef = f$coef[-(1:p)], rows = des$rows))
}

#' Integrated modified OLS (Vogelsang and Wagner, 2014), type 1 regression
#'
#' S_y on S_d, S_x and x. Covariance omega_u.v (Z'Z)^{-1} C'C (Z'Z)^{-1},
#' c_t = S^Z_T - S^Z_{t-1} (partial sums of the regressor matrix Z), with
#' omega_u.v from the FM-OLS long-run quantities of the OLS residuals (as in
#' cointReg). Inference is therefore of the traditional kernel kind, not
#' fixed-b.
#' @noRd
coint_im <- function(y, X, D, kernel, bandwidth) {
  Sy <- cumsum(y)
  Z <- cbind(col_cumsum(D), col_cumsum(X), X)
  f <- ls_fit(Z, Sy)
  S <- col_cumsum(Z)
  ST <- S[nrow(S), ]
  C <- rbind(ST, sweep(-S[-nrow(S), , drop = FALSE], 2, ST, "+"))
  H <- f$ZtZinv %*% crossprod(C) %*% f$ZtZinv
  u_ols <- ls_fit(cbind(D, X), y)$resid
  L <- fm_long_run(u_ols, X, kernel, bandwidth)
  p <- ncol(D) + ncol(X)
  V <- L$omega_uv * H
  list(coef = f$coef[1:p], vcov = V[1:p, 1:p, drop = FALSE], nobs = length(y),
       details = list(omega_uv = L$omega_uv, kernel = kernel,
                      bandwidth = L$lr$bandwidth,
                      bandwidth_rule = L$lr$bandwidth_rule,
                      gamma = f$coef[-(1:p)]))
}

#' Partial sums of every column (a matrix with no columns is returned as is)
#' @noRd
col_cumsum <- function(M) {
  if (ncol(M) == 0) return(M)
  out <- matrix(apply(M, 2, cumsum), nrow = nrow(M))
  colnames(out) <- colnames(M)
  out
}
