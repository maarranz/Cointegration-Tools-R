# Wald tests on the cointegrating vector -------------------------------------------

#' Wald test of linear restrictions on a cointegrating regression
#'
#' Tests \eqn{H_0: R\theta = r} on the coefficients of an FM-OLS, DOLS or
#' IM-OLS fit, with the estimator's own covariance matrix:
#' \deqn{W = (R\hat\theta - r)' (R \hat V R')^{-1} (R\hat\theta - r)
#'   \sim \chi^2(q)}
#' asymptotically. A typical use is \eqn{H_0: \beta = 1}, e.g. a one-for-one
#' long-run relation between two interest rates.
#'
#' Restrictions can be written as text, one equation per element, using the
#' coefficient names of the fit: `"x = 1"`, `"x1 + x2 = 1"`, `"x1 = x2"`.
#' Names that are not syntactic (`(Intercept)`, `log(x)`) go between
#' backticks: ``"`log(x)` = 1"``. Alternatively give `R` and `r` directly.
#'
#' Not available for `method = "OLS"`: its covariance matrix is not valid
#' for a cointegrating regression.
#'
#' @param fit A `"coint_reg"` object.
#' @param hypothesis Character vector of linear restrictions, or a matrix
#'   `R` (one row per restriction, one column per coefficient).
#' @param r Right-hand side when `hypothesis` is a matrix (default zeros).
#' @param alpha Significance level for the decision.
#' @return A `"cointtools_test"` object; `coefficients` holds each
#'   restriction with its estimate \eqn{R\hat\theta}, standard error and
#'   individual z statistic.
#' @examples
#' set.seed(1)
#' n <- 200
#' w <- cumsum(rnorm(n))
#' d <- data.frame(x = w + rnorm(n, sd = 0.4))
#' d$y <- 2 + d$x + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4))
#' fit <- coint_reg(y ~ x, d, method = "FM")
#' coint_wald_test(fit, "x = 1")
#' @export
coint_wald_test <- function(fit, hypothesis, r = NULL, alpha = 0.05) {
  if (!inherits(fit, "coint_reg")) stop("`fit` must be a coint_reg object.", call. = FALSE)
  if (fit$method == "OLS") {
    stop("Wald tests need a covariance that is valid for a cointegrating ",
         "regression: estimate with method = \"FM\", \"D\" or \"IM\".", call. = FALSE)
  }
  b <- fit$coefficients
  V <- fit$vcov
  if (is.character(hypothesis)) {
    rr <- restrictions_from_text(hypothesis, names(b))
    R <- rr$R
    r <- rr$r
    labels <- hypothesis
  } else {
    R <- as.matrix(hypothesis)
    if (ncol(R) != length(b)) {
      stop("`R` must have ", length(b), " columns (one per coefficient: ",
           paste(names(b), collapse = ", "), ").", call. = FALSE)
    }
    r <- r %||% rep(0, nrow(R))
    if (length(r) != nrow(R)) stop("`r` must have one value per row of `R`.", call. = FALSE)
    labels <- restriction_labels(R, r, names(b))
  }
  q <- nrow(R)
  d <- as.numeric(R %*% b - r)
  M <- R %*% V %*% t(R)
  if (qr(M)$rank < q) {
    stop("The restrictions are linearly dependent or not identified.", call. = FALSE)
  }
  w <- as.numeric(t(d) %*% solve(M, d))
  p <- stats::pchisq(w, q, lower.tail = FALSE)
  se <- sqrt(diag(M))
  coefs <- data.frame(restriction = labels, estimate = as.numeric(R %*% b),
                      value = r, std.error = se, z = d / se)
  cov_txt <- switch(fit$method,
    FM = "FM-OLS, omega_u.v (Z'Z)^-1", D = "DOLS, long-run variance of the residuals",
    IM = "IM-OLS, omega_u.v (Z'Z)^-1 C'C (Z'Z)^-1")
  cov_txt <- paste0(cov_txt, "; ", tools::toTitleCase(fit$kernel), " kernel, bandwidth ",
                    format(round(fit$bandwidth, 2), nsmall = 2))
  new_coint_test(list(
    method = "Wald", null_hypothesis = labels, statistic = w,
    distribution = "chi2", df = q, p.value = p, alpha = alpha,
    reject_null = p < alpha, nobs = fit$nobs, coefficients = coefs,
    covariance = cov_txt, R = R, r = r),
    title = paste0("Wald test on the cointegrating regression (",
                   method_labels[[fit$method]], ")"))
}

#' Linear restrictions written as equations, to R and r
#'
#' Each side is a linear expression in the coefficient names. The rows of R
#' are found by evaluating lhs - rhs at the unit vectors; linearity is
#' checked at one more point.
#' @noRd
restrictions_from_text <- function(text, names) {
  k <- length(names)
  rows <- lapply(text, function(h) {
    parts <- strsplit(h, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("Each restriction needs exactly one '=': \"", h, "\".", call. = FALSE)
    }
    ex <- tryCatch(str2lang(paste0("(", parts[1], ") - (", parts[2], ")")),
                   error = function(e) stop("Cannot read the restriction \"", h, "\".",
                                            call. = FALSE))
    used <- all.vars(ex)
    unknown <- setdiff(used, names)
    if (length(unknown)) {
      stop("Unknown coefficient(s) in \"", h, "\": ", paste(unknown, collapse = ", "),
           ". Coefficients: ", paste(names, collapse = ", "),
           " (non-syntactic names between backticks).", call. = FALSE)
    }
    f <- function(v) {
      env <- list2env(stats::setNames(as.list(v), names), parent = baseenv())
      as.numeric(eval(ex, env))
    }
    f0 <- f(rep(0, k))
    Rrow <- vapply(seq_len(k), function(j) f(replace(numeric(k), j, 1)) - f0, numeric(1))
    a <- (-1)^seq_len(k) * (0.37 + 0.71 * seq_len(k))
    if (abs(f(a) - (f0 + sum(Rrow * a))) > 1e-8 * (1 + max(abs(a)))) {
      stop("The restriction \"", h, "\" is not linear in the coefficients.", call. = FALSE)
    }
    if (all(Rrow == 0)) stop("The restriction \"", h, "\" involves no coefficient.", call. = FALSE)
    c(Rrow, -f0)
  })
  M <- do.call(rbind, rows)
  list(R = M[, 1:k, drop = FALSE], r = M[, k + 1])
}

#' Text form of the rows of R theta = r
#' @noRd
restriction_labels <- function(R, r, names) {
  vapply(seq_len(nrow(R)), function(i) {
    nz <- which(R[i, ] != 0)
    terms <- paste0(ifelse(R[i, nz] == 1, "", paste0(format(R[i, nz]), "*")), names[nz])
    paste0(paste(terms, collapse = " + "), " = ", format(r[i]))
  }, character(1))
}
