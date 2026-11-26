# Cointegrating regression: user interface ----------------------------------------

method_labels <- c(OLS = "OLS", FM = "FM-OLS", D = "DOLS", IM = "IM-OLS")
trend_labels <- c(n = "none", c = "constant", ct = "constant and linear trend")

#' Estimate a cointegrating regression
#'
#' Estimates the long-run (levels) regression
#' \deqn{y_t = \alpha + \beta_1 x_{1t} + \dots + \beta_k x_{kt} + z_t}
#' (with a linear trend if `trend = "ct"`) by static OLS, fully modified OLS
#' (FM-OLS), dynamic OLS (DOLS) or integrated modified OLS (IM-OLS). With
#' cointegrated I(1) variables every estimator is superconsistent, but only
#' the last three give standard errors that are valid for inference on
#' \eqn{\beta} when the regressors are endogenous and \eqn{z_t} is serially
#' correlated.
#'
#' * **OLS**: the Engle-Granger step 1 regression. Its conventional standard
#'   errors are not valid here and are not shown; [vcov()] still returns
#'   them, with a warning.
#' * **FM-OLS** (Phillips and Hansen, 1990): OLS corrected for endogeneity
#'   and serial correlation with kernel estimates of the long-run
#'   covariances of the residuals and of \eqn{\Delta x_t}.
#' * **DOLS** (Saikkonen, 1991; Stock and Watson, 1993): OLS with
#'   \eqn{\Delta x_{t-j}}, \eqn{j = -}`leads`, ..., `lags`, added. When `lags`
#'   or `leads` are not given they are chosen by `ic` over 0, ..., `kmax`,
#'   with every pair fitted on the same observations so that the criteria
#'   are comparable; the chosen model is then estimated on its own maximal
#'   sample. If only one of `lags` and `leads` is given, it is held fixed and
#'   the other is selected. Standard errors use the long-run variance of the residuals.
#' * **IM-OLS** (Vogelsang and Wagner, 2014): OLS on the partial sums of the
#'   variables, augmented by \eqn{x_t}; the estimate needs no kernel or
#'   bandwidth. The standard errors use the conditional long-run variance
#'   \eqn{\omega_{u \cdot v}} of FM-OLS (kernel inference, not the fixed-b
#'   inference of the paper).
#'
#' p-values use the standard normal distribution: the results are
#' asymptotic.
#'
#' @param formula A formula `y ~ x1 + x2`, with one variable (or
#'   transformation such as `log(y)`) on each side of every `+`. The
#'   deterministic terms are set by `trend`, not by the formula.
#' @param data A `ts`/`mts`, `zoo`/`xts`, data frame, tibble or tsibble.
#' @param method `"OLS"`, `"FM"`, `"D"` or `"IM"`.
#' @param trend Deterministic terms: `"c"` constant (default), `"ct"`
#'   constant and linear trend, `"n"` none.
#' @param kernel Kernel for the long-run variances: `"bartlett"`,
#'   `"parzen"` or `"qs"` (quadratic spectral).
#' @param bandwidth `"andrews"` (Andrews, 1991, AR(1) rule), `"nw"`
#'   (Newey and West, 1994) or a positive number.
#' @param lags,leads DOLS: numbers of lags and leads of \eqn{\Delta x_t}.
#'   `NULL` selects them by `ic`.
#' @param kmax DOLS: largest number of lags and leads searched; default
#'   \eqn{\lfloor 4 (T/100)^{1/4} \rfloor}.
#' @param ic DOLS: `"aic"` or `"bic"`.
#' @param time Name of the time column of a data frame (detected when there
#'   is exactly one date column).
#' @param start,end Optional first and last observation of the sample, in
#'   the units of the time index (e.g. `c(1990, 1)`, `"1990-01-01"`).
#' @return An object of class `"coint_reg"` with `coefficients`, `vcov`,
#'   `residuals` (the equilibrium error \eqn{\hat z_t} for the whole sample,
#'   as a zoo series), `fitted.values`, `method`, `trend`, `nobs` (the
#'   observations of the estimating regression), `n` (the sample length) and
#'   method details (`omega_uv`, `kernel`, `bandwidth`, `lags`, `leads`,
#'   `selection`, ...).
#' @references
#' Phillips, P. C. B. and Hansen, B. E. (1990). Statistical inference in
#' instrumental variables regression with I(1) processes. *Review of Economic
#' Studies* 57(1), 99-125.
#'
#' Saikkonen, P. (1991). Asymptotically efficient estimation of cointegration
#' regressions. *Econometric Theory* 7(1), 1-21.
#'
#' Stock, J. H. and Watson, M. W. (1993). A simple estimator of cointegrating
#' vectors in higher order integrated systems. *Econometrica* 61(4), 783-820.
#'
#' Vogelsang, T. J. and Wagner, M. (2014). Integrated modified OLS estimation
#' and fixed-b inference for cointegrating regressions. *Journal of
#' Econometrics* 178, 741-760.
#' @seealso [coint_wald_test()]
#' @examples
#' set.seed(1)
#' n <- 200
#' w <- cumsum(rnorm(n))
#' d <- data.frame(x = w + rnorm(n, sd = 0.4))
#' d$y <- 2 + 1.5 * d$x + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4))
#' coint_reg(y ~ x, d, method = "FM")
#' coint_reg(y ~ x, d, method = "D")
#' @export
coint_reg <- function(formula, data, method = c("OLS", "FM", "D", "IM"),
                      trend = c("c", "ct", "n"), kernel = c("bartlett", "parzen", "qs"),
                      bandwidth = "andrews", lags = NULL, leads = NULL,
                      kmax = NULL, ic = c("aic", "bic"), time = NULL,
                      start = NULL, end = NULL) {
  method <- match.arg(method)
  trend <- match.arg(trend)
  kernel <- match.arg(kernel)
  ic <- match.arg(ic)
  check_count(lags, "lags")
  check_count(leads, "leads")
  check_count(kmax, "kmax")

  spec <- parse_coint_formula(formula)
  dat <- as_coint_data(data, data_vars(spec$vars, data, spec$env), time = time,
                       start = start, end = end)
  df <- as.data.frame(zoo::coredata(dat$z))
  y <- eval_term(spec$lhs, df, spec$env)
  X <- do.call(cbind, lapply(spec$rhs, eval_term, df = df, env = spec$env))
  colnames(X) <- spec$rhs_labels
  clash <- intersect(spec$rhs_labels, c("(Intercept)", "trend"))
  if (length(clash)) {
    stop("A regressor is called ", paste0("'", clash, "'", collapse = ", "),
         ", the name of a deterministic term; rename it.", call. = FALSE)
  }
  T <- length(y)
  D <- switch(trend,
              n = matrix(numeric(0), T, 0),
              c = cbind(`(Intercept)` = rep(1, T)),
              ct = cbind(`(Intercept)` = rep(1, T), trend = seq_len(T)))
  if (T < ncol(D) + 2 * ncol(X) + 5) {
    stop("Too few observations (", T, ") for this regression.", call. = FALSE)
  }

  est <- switch(method,
    OLS = coint_ols(y, X, D),
    FM = coint_fm(y, X, D, kernel, bandwidth),
    D = coint_dols(y, X, D, lags, leads, kmax, ic, kernel, bandwidth),
    IM = coint_im(y, X, D, kernel, bandwidth))

  nm <- c(colnames(D), colnames(X))
  coef <- stats::setNames(as.numeric(est$coef), nm)
  V <- est$vcov
  dimnames(V) <- list(nm, nm)
  fit <- as.numeric(cbind(D, X) %*% coef)
  idx <- dat$index

  out <- list(coefficients = coef, vcov = V,
              residuals = zoo::zoo(y - fit, idx),
              fitted.values = zoo::zoo(fit, idx),
              y = zoo::zoo(y, idx), x = zoo::zoo(X, idx),
              method = method, trend = trend, formula = formula,
              yname = spec$lhs_label, xnames = spec$rhs_labels,
              nobs = est$nobs, n = T, index = idx, call = match.call())
  out <- c(out, est$details)
  class(out) <- "coint_reg"
  out
}

#' Check an optional non-negative whole number
#' @noRd
check_count <- function(x, name) {
  if (is.null(x)) return(invisible(TRUE))
  if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x < 0 || x != round(x)) {
    stop("`", name, "` must be a non-negative whole number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Split `y ~ x1 + log(x2)` into its terms
#' @noRd
parse_coint_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be a two-sided formula, e.g. y ~ x1 + x2.", call. = FALSE)
  }
  tt <- stats::terms(formula)
  if (attr(tt, "intercept") == 0) {
    stop("Set the deterministic terms with `trend`, not in the formula ",
         "(use trend = \"n\" for no constant).", call. = FALSE)
  }
  if (any(attr(tt, "order") > 1)) {
    stop("Interactions are not allowed in a cointegrating regression.", call. = FALSE)
  }
  labels <- attr(tt, "term.labels")
  if (length(labels) == 0) stop("The formula has no regressors.", call. = FALSE)
  rhs <- lapply(labels, str2lang)
  env <- environment(formula) %||% parent.frame()
  list(lhs = formula[[2]], lhs_label = deparse1(formula[[2]]),
       rhs = rhs, rhs_labels = labels, vars = all.vars(formula), env = env)
}

#' Evaluate one term on the data, as a numeric vector
#' @noRd
eval_term <- function(term, df, env = environment(stats::lm)) {
  v <- suppressWarnings(eval(term, df, env))
  if (!is.numeric(v) || length(v) != nrow(df) || !is.null(dim(v)) && ncol(v) != 1) {
    stop("The term `", deparse1(term), "` does not give one numeric series.", call. = FALSE)
  }
  if (anyNA(v) || any(!is.finite(v))) {
    stop("The term `", deparse1(term), "` has missing or infinite values ",
         "(e.g. the log of a non-positive number).", call. = FALSE)
  }
  as.numeric(v)
}

#' @export
coef.coint_reg <- function(object, ...) object$coefficients

#' @export
vcov.coint_reg <- function(object, ...) {
  if (object$method == "OLS") {
    warning("Conventional OLS covariance: not valid for inference on a ",
            "cointegrating vector; use method = \"FM\", \"D\" or \"IM\".", call. = FALSE)
  }
  object$vcov
}

#' @export
residuals.coint_reg <- function(object, ...) object$residuals

#' @export
fitted.coint_reg <- function(object, ...) object$fitted.values

#' @export
nobs.coint_reg <- function(object, ...) object$nobs

#' Coefficient table of a cointegrating regression
#' @noRd
coint_table <- function(object, level = 0.95) {
  b <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- b / se
  q <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(term = names(b), estimate = unname(b), std.error = unname(se),
             statistic = unname(z), p.value = unname(2 * stats::pnorm(-abs(z))),
             conf.low = unname(b - q * se), conf.high = unname(b + q * se),
             row.names = NULL)
}

#' Summary of a cointegrating regression
#'
#' @param object A `"coint_reg"` object.
#' @param level Confidence level of the intervals.
#' @param ... Not used.
#' @return An object of class `"summary.coint_reg"`: the fit with a
#'   coefficient table (estimates, standard errors, z statistics, normal
#'   p-values and confidence intervals). For OLS the table has the estimates
#'   only.
#' @export
summary.coint_reg <- function(object, level = 0.95, ...) {
  check_prob(level, "level")
  tab <- coint_table(object, level)
  if (object$method == "OLS") tab <- tab[, c("term", "estimate")]
  structure(c(unclass(object), list(table = tab, level = level)),
            class = "summary.coint_reg")
}

#' @export
print.coint_reg <- function(x, digits = 4, ...) {
  print(summary(x), digits = digits)
  invisible(x)
}

#' @export
print.summary.coint_reg <- function(x, digits = 4, ...) {
  title <- paste0("Cointegrating regression (", method_labels[[x$method]], ")")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Model:         ", x$yname, " ~ ", paste(x$xnames, collapse = " + "), "\n", sep = "")
  cat("Deterministic: ", trend_labels[[x$trend]], "\n", sep = "")
  cat("Sample:        ", format(x$index[1]), " to ", format(x$index[x$n]),
      " (T = ", x$n, if (x$nobs != x$n) paste0("; ", x$nobs, " in the regression"),
      ")\n", sep = "")
  if (x$method == "D") {
    sel <- x$selection
    how <- if (is.null(sel)) "given" else {
      free <- c("lags", "leads")[!sel$fixed]
      paste0(paste(free, collapse = " and "), " by ", toupper(sel$ic), " over 0..", sel$kmax,
             if (any(sel$fixed)) paste0(", ", c("lags", "leads")[sel$fixed], " fixed"),
             ", on a common sample of ", sel$nobs, " obs.")
    }
    cat("Leads, lags:   ", x$leads, " leads, ", x$lags, " lags of the differenced regressors (", how,
        ")\n", sep = "")
  }
  if (x$method != "OLS") {
    rule <- switch(x$bandwidth_rule, andrews = "Andrews", nw = "Newey-West", fixed = "fixed")
    lrv <- if (x$method == "D") "long-run variance of the residuals"
           else "conditional long-run variance omega_u.v"
    cat("Inference:     ", lrv, ", ", tools::toTitleCase(x$kernel), " kernel, bandwidth ",
        format(round(x$bandwidth, 2), nsmall = 2), " (", rule, ")\n", sep = "")
  }
  cat("\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  if (x$method != "OLS") {
    names(tab) <- c("", "Estimate", "Std. Error", "z", "P(>|z|)",
                    paste0(format(100 * (1 - x$level) / 2), "%"),
                    paste0(format(100 * (1 + x$level) / 2), "%"))
  } else {
    names(tab) <- c("", "Estimate")
  }
  print(tab, row.names = FALSE)
  if (x$method == "OLS") {
    cat("\nOLS standard errors are not valid for a cointegrating regression;",
        "\nuse method = \"FM\", \"D\" or \"IM\" for inference on the coefficients.\n")
  }
  invisible(x)
}

#' Check a probability-type argument (significance or confidence level)
#' @noRd
check_prob <- function(x, name) {
  if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", name, "` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Variables of a formula to take from the data: all, except those that are
#' not columns of the data but exist in the formula's environment (constants
#' or objects used by a transformation)
#' @noRd
data_vars <- function(vars, data, env) {
  cols <- colnames(data) %||% names(data)
  keep <- vars %in% cols | !vapply(vars, exists, logical(1), envir = env)
  vars[keep]
}
