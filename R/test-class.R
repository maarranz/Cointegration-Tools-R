# Test results ----------------------------------------------------------------
#
# Every test returns a list of class "cointtools_test" with the fields of the
# suite (adltools_test, forecastingtools_test): method, null_hypothesis,
# statistic, distribution ("normal", "chi2", ...), df, p.value, alpha,
# reject_null, nobs, coefficients (data frame or NULL); plus `covariance`, a
# one-line description of the covariance used. Tests can add their own fields.

#' Build a test result
#' @noRd
new_coint_test <- function(x, title) {
  defaults <- list(alternative = NULL, covariance = NULL, coefficients = NULL,
                   df = NA_integer_, critical = NULL)
  for (nm in names(defaults)) if (!nm %in% names(x)) x[nm] <- list(defaults[[nm]])
  structure(x, class = "cointtools_test", title = title)
}

fmt_num <- function(x, digits) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("-")
  if (is.numeric(x) && !is.integer(x)) format(signif(x, digits)) else format(x)
}

#' Print a test result
#'
#' Prints the hypothesis, statistic, p-value and decision, followed by the
#' coefficient table when there is one. Every element can also be extracted
#' directly, e.g. `result$p.value`.
#'
#' @param x A test result from cointegrationtools.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.cointtools_test <- function(x, digits = 4, ...) {
  title <- attr(x, "title")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  label <- switch(x$distribution, normal = "z", chi2 = "Chi-square", x$distribution)
  df <- if (!is.null(x$df) && !is.na(x$df)) paste0("(", x$df, ")") else ""
  cat("H0:          ", paste(x$null_hypothesis, collapse = "; "), "\n", sep = "")
  if (!is.null(x$alternative)) cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("Statistic:   ", fmt_num(x$statistic, digits), "  [", label, df, "]\n", sep = "")
  if (!is.null(x$critical)) {
    cat("Critical:    ", fmt_num(x$critical, digits), " at ", format(100 * x$alpha),
        "%\n", sep = "")
  }
  cat("p-value:     ", fmt_num(x$p.value, digits), "\n", sep = "")
  cat("Decision:    ", if (isTRUE(x$reject_null)) "reject H0" else "do not reject H0",
      " at ", format(100 * x$alpha), "%\n", sep = "")
  if (!is.null(x$covariance)) cat("Covariance:  ", x$covariance, "\n", sep = "")
  cat("N:           ", x$nobs, "\n", sep = "")
  if (is.data.frame(x$coefficients) && nrow(x$coefficients)) {
    cat("\n")
    tab <- x$coefficients
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}
