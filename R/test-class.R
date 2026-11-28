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
                   df = NA_integer_, critical_values = NULL, details = NULL,
                   warnings = character(), p_bounded = FALSE, assumption = NULL)
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
  if (length(x$details)) {
    for (i in seq_along(x$details)) {
      cat(formatC(paste0(names(x$details)[i], ":"), width = -13), x$details[[i]], "\n", sep = "")
    }
    cat("\n")
  }
  cat("H0:          ", paste(x$null_hypothesis, collapse = "; "), "\n", sep = "")
  if (!is.null(x$alternative)) cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("Statistic:   ", fmt_num(x$statistic, digits), "  [", label, df, "]\n", sep = "")
  if (length(x$critical_values)) {
    cv <- x$critical_values
    cat("Critical:    ", paste0(format(100 * as.numeric(names(cv))), "% ",
                                format(round(cv, 3), nsmall = 3), collapse = "   "), "\n", sep = "")
  }
  # the bounds of the simulated tables apply only to simulated p-values (L01)
  pv <- if (!is.null(x$p_text)) x$p_text
        else if (isTRUE(x$p_bounded) && x$p.value < 0.5) "< 0.001"
        else if (isTRUE(x$p_bounded)) "> 0.999" else fmt_num(x$p.value, digits)
  cat("p-value:     ", pv, "\n", sep = "")
  dec <- if (is.na(x$reject_null)) "undetermined (p is only known to be below 0.001)"
         else if (x$reject_null) "reject H0" else "do not reject H0"
  cat("Decision:    ", dec, " at ", format(100 * x$alpha), "%\n", sep = "")
  if (!is.null(x$assumption)) cat("Assumes:     ", x$assumption, "\n", sep = "")
  if (!is.null(x$covariance)) cat("Covariance:  ", x$covariance, "\n", sep = "")
  cat("N:           ", x$nobs, "\n", sep = "")
  if (is.data.frame(x$statistics)) {
    cat("\nAll three statistics (minimum over the break dates):\n")
    st <- x$statistics
    st$value <- signif(st$value, digits)
    print(st, row.names = FALSE)
  }
  if (is.data.frame(x$normalisations)) {
    cat("\nEach variable as the dependent one:\n")
    nt <- x$normalisations
    nt$statistic <- signif(nt$statistic, digits)
    nt$p.value <- signif(nt$p.value, digits)
    print(nt, row.names = FALSE)
  }
  for (w in x$warnings) cat("Note: ", w, "\n", sep = "")
  if (is.data.frame(x$coefficients) && nrow(x$coefficients)) {
    cat("\n")
    tab <- x$coefficients
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}
