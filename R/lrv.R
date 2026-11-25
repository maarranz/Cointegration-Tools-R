# Long-run covariance matrices ------------------------------------------------
#
# Kernel estimators of the long-run covariance Omega and of the one-sided
# long-run covariance Delta of a stationary vector series w_t (T x m):
#   Gamma_j = (1/T) sum_{t>j} w_t w_{t-j}'           (j >= 0)
#   Omega   = Gamma_0 + sum_{j>=1} k(j/b) (Gamma_j + Gamma_j')
#   Delta   = Gamma_0 + sum_{j>=1} k(j/b) Gamma_j'    (= sum_{j>=0} E w_{t-j} w_t')
# Kernels: Bartlett, Parzen, quadratic spectral. Bandwidths: Andrews (1991)
# with AR(1) approximations for every column (equal weights), Newey and West
# (1994), or a number. The conventions follow cointReg (Aschersleben and
# Wagner, 2016), so that the estimators can be checked against it, except
# that a Bartlett or Parzen bandwidth below 1 gives no autocovariances (cointReg
# then uses lag 1 with a negative weight).

kernel_names <- c("bartlett", "parzen", "qs")

#' Kernel weights k(j / b) for j = 1, ..., T - 1
#' @noRd
kernel_weights <- function(n, b, kernel) {
  j <- seq_len(n - 1)
  x <- j / b
  switch(kernel,
    bartlett = ifelse(x < 1, 1 - x, 0),
    parzen = ifelse(x <= 0.5, 1 - 6 * x^2 + 6 * x^3,
                    ifelse(x < 1, 2 * (1 - x)^3, 0)),
    qs = {
      a <- 6 * pi * x / 5
      25 / (12 * pi^2 * x^2) * (sin(a) / a - cos(a))
    })
}

#' Andrews (1991) automatic bandwidth, AR(1) approximation for each column
#' @noRd
andrews_bandwidth <- function(w, kernel) {
  w <- as.matrix(w)
  n <- nrow(w)
  rho <- s2 <- numeric(ncol(w))
  for (i in seq_len(ncol(w))) {
    u0 <- w[-1, i]
    u1 <- w[-n, i]
    rho[i] <- sum(u0 * u1) / sum(u1^2)
    s2[i] <- sum((u0 - rho[i] * u1)^2) / n
  }
  den <- sum(s2^2 / (1 - rho)^4)
  if (kernel == "bartlett") {
    a <- sum(4 * rho^2 * s2^2 / ((1 - rho)^6 * (1 + rho)^2)) / den
    b <- 1.1447 * (a * n)^(1 / 3)
  } else {
    a <- sum(4 * rho^2 * s2^2 / (1 - rho)^8) / den
    b <- switch(kernel, parzen = 2.6614, qs = 1.3221) * (a * n)^(1 / 5)
  }
  min(b, n - 1)
}

#' Newey and West (1994) automatic bandwidth (columns summed, equal weights)
#' @noRd
newey_west_bandwidth <- function(w, kernel) {
  w <- as.matrix(w)
  n <- nrow(w)
  pw <- switch(kernel, bartlett = 2 / 9, parzen = 4 / 25, qs = 2 / 25)
  m <- floor(4 * (n / 100)^pw)
  u <- rowSums(w)
  s <- vapply(0:m, function(j) sum(u[(j + 1):n] * u[1:(n - j)]) / n, numeric(1))
  s0 <- s[1] + 2 * sum(s[-1])
  if (kernel == "bartlett") {
    q <- 1
    sq <- 2 * sum(seq_len(m) * s[-1])
  } else {
    q <- 2
    sq <- 2 * sum(seq_len(m)^2 * s[-1])
  }
  cst <- switch(kernel, bartlett = 1.1447, parzen = 2.6614, qs = 1.3221)
  cst * ((sq / s0)^2)^(1 / (2 * q + 1)) * n^(1 / (2 * q + 1))
}

#' Resolve the bandwidth argument
#' @noRd
resolve_bandwidth <- function(bandwidth, w, kernel) {
  if (is.numeric(bandwidth)) {
    if (length(bandwidth) != 1 || !is.finite(bandwidth) || bandwidth <= 0) {
      stop("A numeric `bandwidth` must be a single positive number.", call. = FALSE)
    }
    return(list(value = bandwidth, rule = "fixed"))
  }
  rule <- match.arg(bandwidth, c("andrews", "nw"))
  value <- if (rule == "andrews") andrews_bandwidth(w, kernel)
           else newey_west_bandwidth(w, kernel)
  list(value = value, rule = rule)
}

#' Long-run covariance (Omega), one-sided (Delta) and Gamma_0 of w
#' @noRd
long_run_cov <- function(w, kernel = "bartlett", bandwidth = "andrews") {
  w <- as.matrix(w)
  n <- nrow(w)
  bw <- resolve_bandwidth(bandwidth, w, kernel)
  k <- kernel_weights(n, bw$value, kernel)
  G0 <- crossprod(w) / n
  Omega <- Delta <- G0
  for (j in which(k != 0)) {
    Gj <- crossprod(w[(j + 1):n, , drop = FALSE], w[1:(n - j), , drop = FALSE]) / n
    Omega <- Omega + k[j] * (Gj + t(Gj))
    Delta <- Delta + k[j] * t(Gj)
  }
  dimnames(Omega) <- dimnames(Delta) <- dimnames(G0) <- list(colnames(w), colnames(w))
  list(Omega = Omega, Delta = Delta, Sigma = G0, bandwidth = bw$value,
       bandwidth_rule = bw$rule, kernel = kernel)
}

#' Split a long-run matrix of (u, v) into its blocks (u is the first column)
#' @noRd
lr_blocks <- function(M) {
  m <- ncol(M)
  list(uu = M[1, 1], uv = M[1, 2:m, drop = FALSE], vu = M[2:m, 1, drop = FALSE],
       vv = M[2:m, 2:m, drop = FALSE])
}
