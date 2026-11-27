# Recursive and rolling estimation of the cointegrating regression --------------

#' Recursive or rolling estimates of a cointegrating regression
#'
#' Re-estimates the cointegrating regression of [coint_reg()] on a growing
#' sample (recursive: observations 1, ..., t for t = `min_obs`, ..., T) or on
#' a moving window of fixed length (rolling: t - `window` + 1, ..., t), and
#' collects the coefficients with their confidence intervals. Plotting them
#' against t shows whether the long-run relation is stable: a coefficient
#' that drifts, or jumps, outside the bands of earlier estimates points to a
#' break (and, for a Gregory-Hansen test later, to its date).
#'
#' Each estimate uses only its own sample, including the automatic bandwidth
#' and, with `trend = "ct"`, a trend that starts at 1 in every window. For
#' DOLS, leads and lags not given are chosen once on the whole sample and
#' then held fixed, so that the estimates differ only through the data.
#' For OLS no intervals are reported (its standard errors are not valid).
#'
#' @inheritParams coint_reg
#' @param window `NULL` (default) for recursive estimation, or the length of
#'   a rolling window.
#' @param min_obs Smallest sample for recursive estimation (default: the
#'   larger of 30 and 15% of the sample).
#' @param level Confidence level of the intervals.
#' @return A data frame of class `"coint_recursive"` in long form: `end`
#'   (the last observation of each sample, in the units of the time index),
#'   `nobs`, `term`, `estimate`, `std.error`, `conf.low`, `conf.high`; the
#'   whole-sample fit is in `attr(x, "full")`. `plot()` draws one panel per
#'   coefficient.
#' @seealso [coint_reg()], [plot.coint_recursive()]
#' @examples
#' rec <- coint_recursive(gs10 ~ tb3ms, term_structure, method = "FM")
#' plot(rec, ref = c(tb3ms = 1))
#' roll <- coint_recursive(gs10 ~ tb3ms, term_structure, method = "FM", window = 120)
#' plot(roll, terms = "tb3ms")
#' @export
coint_recursive <- function(formula, data, method = c("FM", "D", "IM", "OLS"),
                            trend = c("c", "ct", "n"), window = NULL, min_obs = NULL,
                            kernel = c("bartlett", "parzen", "qs"), bandwidth = "andrews",
                            lags = NULL, leads = NULL, kmax = NULL, ic = c("aic", "bic"),
                            level = 0.95, time = NULL, start = NULL, end = NULL) {
  method <- match.arg(method)
  trend <- match.arg(trend)
  kernel <- match.arg(kernel)
  ic <- match.arg(ic)
  check_prob(level, "level")
  check_count(window, "window")
  check_count(min_obs, "min_obs")
  full <- coint_reg(formula, data, method = method, trend = trend, kernel = kernel,
                    bandwidth = bandwidth, lags = lags, leads = leads, kmax = kmax,
                    ic = ic, time = time, start = start, end = end)
  if (method == "D") {
    lags <- full$lags
    leads <- full$leads
  }
  y <- as.numeric(zoo::coredata(full$y))
  X <- as.matrix(zoo::coredata(full$x))
  colnames(X) <- full$xnames
  T <- full$n
  k <- ncol(X)
  m <- c(n = 0, c = 1, ct = 2)[[trend]]
  smallest <- m + 2 * k + 10 + if (method == "D") k * (lags + leads + 2) + lags + leads else 0
  if (is.null(window)) {
    min_obs <- min_obs %||% max(30L, ceiling(0.15 * T), smallest)
    if (min_obs >= T) stop("`min_obs` must be smaller than the sample (", T, ").", call. = FALSE)
    ends <- seq.int(min_obs, T)
    first <- rep(1L, length(ends))
  } else {
    if (window < smallest || window >= T) {
      stop("`window` must be between ", smallest, " and ", T - 1, ".", call. = FALSE)
    }
    ends <- seq.int(window, T)
    first <- ends - window + 1L
  }
  if (min(ends - first + 1L) < smallest) {
    stop("Some samples are too short for this regression (fewer than ", smallest,
         " observations).", call. = FALSE)
  }
  q <- stats::qnorm(1 - (1 - level) / 2)
  nm <- names(full$coefficients)
  rows <- lapply(seq_along(ends), function(i) {
    r <- first[i]:ends[i]
    n <- length(r)
    D <- switch(trend, n = matrix(numeric(0), n, 0), c = matrix(1, n, 1),
                ct = cbind(1, seq_len(n)))
    est <- switch(method,
      OLS = coint_ols(y[r], X[r, , drop = FALSE], D),
      FM = coint_fm(y[r], X[r, , drop = FALSE], D, kernel, bandwidth),
      D = coint_dols(y[r], X[r, , drop = FALSE], D, lags, leads, NULL, ic, kernel, bandwidth),
      IM = coint_im(y[r], X[r, , drop = FALSE], D, kernel, bandwidth))
    b <- as.numeric(est$coef)
    se <- if (method == "OLS") rep(NA_real_, length(b)) else sqrt(diag(est$vcov))
    data.frame(i = ends[i], nobs = est$nobs, term = nm, estimate = b, std.error = se,
               conf.low = b - q * se, conf.high = b + q * se)
  })
  out <- do.call(rbind, rows)
  out$end <- full$index[out$i]
  out <- out[c("end", "nobs", "term", "estimate", "std.error", "conf.low", "conf.high")]
  rownames(out) <- NULL
  structure(out, class = c("coint_recursive", "data.frame"), full = full,
            type = if (is.null(window)) "recursive" else "rolling", window = window,
            level = level, method = method)
}

#' @export
print.coint_recursive <- function(x, digits = 4, ...) {
  full <- attr(x, "full")
  type <- attr(x, "type")
  title <- paste0(tools::toTitleCase(type), " estimates of the cointegrating regression (",
                  method_labels[[attr(x, "method")]], ")")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Model:   ", full$yname, " ~ ", paste(full$xnames, collapse = " + "),
      " (", trend_labels[[full$trend]], ")\n", sep = "")
  ends <- unique(x$end)
  cat("Samples: ", length(ends), if (type == "rolling") paste0(" windows of ", attr(x, "window"), " obs.")
      else " growing samples", ", last observation from ", format(ends[1]), " to ",
      format(ends[length(ends)]), "\n\n", sep = "")
  tab <- do.call(rbind, lapply(split(x, factor(x$term, levels = unique(x$term))), function(d) {
    data.frame(term = d$term[1], first = d$estimate[1], min = min(d$estimate),
               max = max(d$estimate), last = d$estimate[nrow(d)])
  }))
  tab[-1] <- lapply(tab[-1], signif, digits = digits)
  print(tab, row.names = FALSE)
  cat("\nUse plot() to see the estimates through time.\n")
  invisible(x)
}

#' Plot recursive or rolling estimates
#'
#' One panel per coefficient: the estimates against the last observation of
#' each sample, with their confidence band (not for OLS), the whole-sample
#' estimate (dashed) and, optionally, reference values such as a
#' theoretical coefficient.
#'
#' @param x A `"coint_recursive"` object.
#' @param terms Coefficients to plot (default: all).
#' @param ref Optional named numeric vector of reference values, e.g.
#'   `c(tb3ms = 1)`, drawn as dotted lines.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.coint_recursive <- function(x, terms = NULL, ref = NULL, ...) {
  terms <- terms %||% unique(x$term)
  bad <- setdiff(terms, x$term)
  if (length(bad)) stop("Unknown term(s): ", paste(bad, collapse = ", "), call. = FALSE)
  full <- attr(x, "full")
  op <- graphics::par(mfrow = c(length(terms), 1), mar = c(3, 4, 2.5, 1))
  on.exit(graphics::par(op))
  for (tm in terms) {
    d <- x[x$term == tm, ]
    rng <- range(c(d$estimate, d$conf.low, d$conf.high, ref[tm]), na.rm = TRUE)
    graphics::plot(d$end, d$estimate, type = "n", ylim = rng, xlab = "", ylab = tm,
                   main = paste0(tm, ": ", attr(x, "type"), " ", method_labels[[attr(x, "method")]],
                                 " estimate"), ...)
    if (!all(is.na(d$conf.low))) {
      graphics::polygon(c(d$end, rev(d$end)), c(d$conf.low, rev(d$conf.high)),
                        col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
    }
    graphics::lines(d$end, d$estimate, lwd = 2, col = "steelblue4")
    graphics::abline(h = full$coefficients[[tm]], lty = 2)
    if (!is.null(ref) && tm %in% names(ref)) graphics::abline(h = ref[[tm]], lty = 3, col = "red")
  }
  invisible(x)
}
