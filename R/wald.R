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
#' backticks: ``"`log(x)` = 1"``. Only linear expressions are accepted:
#' numbers, coefficient names, `+`, `-`, parentheses, and products or
#' quotients by numbers (`"2*x1 - x2/3 = 0"`). Alternatively give `R` and
#' `r` directly; if `R` has column names they are matched to the coefficient
#' names, otherwise columns are taken in the order of `coef(fit)`.
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
  check_prob(alpha, "alpha")
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
    R <- hypothesis
    if (is.null(dim(R))) R <- matrix(R, nrow = 1, dimnames = list(NULL, names(R)))
    R <- as.matrix(R)
    if (!is.numeric(R) || ncol(R) != length(b)) {
      stop("`R` must be numeric with ", length(b), " columns (one per coefficient: ",
           paste(names(b), collapse = ", "), ").", call. = FALSE)
    }
    if (!is.null(colnames(R))) {
      # named columns are matched to the coefficients by name, not position
      if (!setequal(colnames(R), names(b)) || anyDuplicated(colnames(R))) {
        stop("The column names of `R` must be the coefficient names: ",
             paste(names(b), collapse = ", "), ".", call. = FALSE)
      }
      R <- R[, names(b), drop = FALSE]
    }
    r <- r %||% rep(0, nrow(R))
    if (length(r) != nrow(R)) stop("`r` must have one value per row of `R`.", call. = FALSE)
    if (any(!is.finite(R)) || any(!is.finite(r))) {
      stop("`R` and `r` must be finite.", call. = FALSE)
    }
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
#' Each side is parsed (never evaluated) by linear_form(): only numbers,
#' coefficient names, parentheses, + and -, and products or quotients by
#' numbers are accepted, so a nonlinear restriction or any other R call is
#' refused.
#' @noRd
restrictions_from_text <- function(text, names) {
  k <- length(names)
  rows <- lapply(text, function(h) {
    parts <- strsplit(h, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2 || !nzchar(trimws(parts[1])) || !nzchar(trimws(parts[2]))) {
      stop("Each restriction needs exactly one '=': \"", h, "\".", call. = FALSE)
    }
    side <- function(txt) {
      ex <- tryCatch(str2lang(txt), error = function(e)
        stop("Cannot read the restriction \"", h, "\".", call. = FALSE))
      linear_form(ex, names, h)
    }
    l <- side(parts[1])
    r <- side(parts[2])
    Rrow <- l$coef - r$coef
    if (all(Rrow == 0)) stop("The restriction \"", h, "\" involves no coefficient.", call. = FALSE)
    c(Rrow, r$const - l$const)
  })
  M <- do.call(rbind, rows)
  list(R = M[, 1:k, drop = FALSE], r = M[, k + 1])
}

#' A linear expression as coefficients on `names` plus a constant
#' @noRd
linear_form <- function(ex, names, h) {
  k <- length(names)
  bad <- function(why) stop("The restriction \"", h, "\" is not a linear restriction on ",
                            "the coefficients (", why, "). Coefficients: ",
                            paste(names, collapse = ", "),
                            " (non-syntactic names between backticks).", call. = FALSE)
  lin <- function(coef = numeric(k), const = 0) list(coef = coef, const = const)
  rec <- function(e) {
    if (is.numeric(e) && length(e) == 1 && is.finite(e)) return(lin(const = e))
    if (is.symbol(e)) {
      nm <- as.character(e)
      j <- match(nm, names)
      if (is.na(j)) bad(paste0("unknown coefficient '", nm, "'"))
      return(lin(replace(numeric(k), j, 1)))
    }
    if (!is.call(e)) bad("unsupported element")
    op <- as.character(e[[1]])
    args <- as.list(e)[-1]
    if (op == "(" && length(args) == 1) return(rec(args[[1]]))
    if (op %in% c("+", "-") && length(args) == 1) {
      a <- rec(args[[1]])
      return(if (op == "-") lin(-a$coef, -a$const) else a)
    }
    if (op %in% c("+", "-") && length(args) == 2) {
      a <- rec(args[[1]]); b <- rec(args[[2]])
      sg <- if (op == "-") -1 else 1
      return(lin(a$coef + sg * b$coef, a$const + sg * b$const))
    }
    if (op == "*" && length(args) == 2) {
      a <- rec(args[[1]]); b <- rec(args[[2]])
      if (any(a$coef != 0) && any(b$coef != 0)) bad("a product of coefficients")
      if (any(a$coef != 0)) return(lin(a$coef * b$const, a$const * b$const))
      return(lin(b$coef * a$const, a$const * b$const))
    }
    if (op == "/" && length(args) == 2) {
      a <- rec(args[[1]]); b <- rec(args[[2]])
      if (any(b$coef != 0)) bad("division by a coefficient")
      if (b$const == 0) bad("division by zero")
      return(lin(a$coef / b$const, a$const / b$const))
    }
    bad(paste0("'", op, "' is not allowed"))
  }
  rec(ex)
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
