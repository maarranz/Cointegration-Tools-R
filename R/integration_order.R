# Order of integration of each variable -------------------------------------------------

#' Order of integration of the variables of a cointegrating regression
#'
#' Step 0 of a cointegration analysis: cointegration is only a question for
#' variables that are I(1). For each variable, a unit-root specification
#' strategy (`unitroottools::ur_strategy()`) is applied to the first
#' difference and to the level, and read from the highest order down
#' (Dickey and Pantula, 1987):
#'
#' * unit root not rejected in the difference: **I(2) or higher** (the case
#'   of the price-level example in the chapter), whatever the level test
#'   says: testing one unit root in an I(2) series is not a valid test;
#' * rejected in the difference, not in the level: **I(1)**;
#' * rejected in both: **I(0)**.
#'
#' Lags are chosen by BIC by default. MAIC (the default of
#' `ur_strategy()`) tends to choose very long lags for series that are far
#' from a unit root, such as the differences of an I(1) series, and then
#' misses the rejection.
#'
#' With `strategy = "hamilton"`, `trending` refers to the levels; the
#' differences are always tested as non-trending series (the difference of a
#' series with drift has a constant mean).
#'
#' @param formula A formula naming the variables, with or without a
#'   left-hand side: `y ~ x1 + x2` (as in [eg_test()]) or `~ y + x1 + x2`.
#'   Transformations such as `log(y)` are allowed.
#' @inheritParams coint_reg
#' @param strategy,trending,lags,max_lags,criterion Passed to
#'   `unitroottools::ur_strategy()`. `trending` can be a single value or one
#'   per variable (in the order of the formula).
#' @param level Significance level of the tests.
#' @return A data frame of class `"integration_order"`, one row per variable:
#'   the specification chosen, statistic and p-value for the level and for
#'   the difference, and the conclusion (`order`). The full strategy results
#'   are in `attr(x, "strategies")`.
#' @references
#' Dickey, D. A. and Pantula, S. G. (1987). Determining the order of
#' differencing in autoregressive processes. *Journal of Business & Economic
#' Statistics* 5(4), 455-461.
#' @seealso `unitroottools::ur_strategy()`, [eg_test()]
#' @examples
#' set.seed(1)
#' n <- 200
#' w <- cumsum(rnorm(n))
#' d <- data.frame(x = w + rnorm(n, sd = 0.4), u = rnorm(n))
#' d$y <- 2 + 1.5 * d$x + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4))
#' integration_order(~ y + x + u, d)
#' @export
integration_order <- function(formula, data, strategy = c("original", "conservative", "hamilton"),
                              trending = NULL, lags = NULL, max_lags = NULL,
                              criterion = c("BIC", "AIC", "MAIC", "MBIC"), level = 0.05,
                              time = NULL, start = NULL, end = NULL) {
  strategy <- match.arg(strategy)
  criterion <- match.arg(criterion)
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula, e.g. ~ y + x or y ~ x.", call. = FALSE)
  }
  terms <- c(if (length(formula) == 3) deparse1(formula[[2]]),
             attr(stats::terms(formula), "term.labels"))
  if (!length(terms)) stop("The formula names no variable.", call. = FALSE)
  if (!is.null(trending) && !length(trending) %in% c(1L, length(terms))) {
    stop("`trending` must be a single value or one per variable (", length(terms), ").",
         call. = FALSE)
  }
  trending <- if (is.null(trending)) vector("list", length(terms)) else
    as.list(rep_len(trending, length(terms)))

  dat <- as_coint_data(data, all.vars(formula), time = time, start = start, end = end)
  df <- as.data.frame(zoo::coredata(dat$z))
  run <- function(y, tr) {
    unitroottools::ur_strategy(y, strategy = strategy, trending = tr, lags = lags,
                               max_lags = max_lags, criterion = criterion, level = level)
  }
  strategies <- list()
  rows <- lapply(seq_along(terms), function(i) {
    y <- eval_term(str2lang(terms[i]), df)
    lev <- run(y, trending[[i]])
    dif <- run(diff(y), if (strategy == "hamilton") FALSE else NULL)
    strategies[[terms[i]]] <<- list(level = lev, difference = dif)
    order <- if (!dif$unit_root_rejected) "I(2) or higher"
             else if (lev$unit_root_rejected) "I(0)" else "I(1)"
    data.frame(variable = terms[i],
               level_spec = lev$specification, level_stat = lev$terminal$statistic,
               level_p = lev$terminal$p_value,
               diff_spec = dif$specification, diff_stat = dif$terminal$statistic,
               diff_p = dif$terminal$p_value, order = order)
  })
  out <- do.call(rbind, rows)
  structure(out, class = c("integration_order", "data.frame"), strategies = strategies,
            strategy = strategy, level = level, n = nrow(df))
}

#' @export
print.integration_order <- function(x, digits = 3, ...) {
  title <- "Order of integration"
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  cat("Unit-root strategy: ", attr(x, "strategy"), " (unitroottools), level ",
      format(100 * attr(x, "level")), "%; T = ", attr(x, "n"), "\n\n", sep = "")
  tab <- data.frame(Variable = x$variable,
                    `Level: spec` = x$level_spec,
                    `stat` = round(x$level_stat, digits),
                    `p` = fmt_p(x$level_p),
                    `Difference: spec` = x$diff_spec,
                    `stat ` = round(x$diff_stat, digits),
                    `p ` = fmt_p(x$diff_p),
                    Order = x$order, check.names = FALSE)
  print(tab, row.names = FALSE)
  cat("\nspec: deterministic terms chosen by the strategy (n none, c constant,",
      "ct constant and trend).\n")
  if (any(x$order != "I(1)")) {
    cat("Cointegration tests assume that every variable is I(1).\n")
  }
  invisible(x)
}

#' p-values with the bounds of the simulated tables
#' @noRd
fmt_p <- function(p) {
  ifelse(is.na(p), "-", ifelse(p <= 0.001 + 1e-9, "<0.001",
         ifelse(p >= 0.999 - 1e-9, ">0.999", formatC(p, format = "f", digits = 3))))
}
