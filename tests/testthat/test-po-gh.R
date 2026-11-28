test_that("po_stats follows the Phillips-Perron formulas", {
  set.seed(41)
  u <- cumsum(rnorm(150)) * 0.2 + arima.sim(list(ar = 0.6), 150)
  u <- as.numeric(u)
  ps <- cointegrationtools:::po_stats(u, "bartlett", 4)
  T <- 150; N <- T - 1
  r <- lm(u[-1] ~ 0 + u[-T])
  rho <- coef(r)[[1]]; e <- residuals(r)
  s2 <- sum(e^2) / N
  g <- sapply(1:3, function(j) sum(e[(j + 1):N] * e[1:(N - j)]) / N)
  l2 <- s2 + 2 * sum((1 - (1:3) / 4) * g)
  S <- sum(u[-T]^2) / N^2
  tr <- coef(summary(r))[1, 1] - 1
  tr <- (rho - 1) / coef(summary(r))[1, 2]
  expect_equal(ps$Za, N * (rho - 1) - (l2 - s2) / (2 * S))
  expect_equal(ps$Zt, sqrt(s2 / l2) * tr - (l2 - s2) / (2 * sqrt(l2 * S)))
})

test_that("without autocovariances Z_t is the Engle-Granger statistic with no lags", {
  d <- simulate_coint(120, rho = 0.6, seed = 42)
  for (tr in c("n", "c", "ct")) {
    p <- po_test(y ~ x, d, trend = tr, bandwidth = 0.5)
    e <- suppressMessages(eg_test(y ~ x, d, trend = tr, lags = 0))
    expect_equal(p$statistics[["Zt"]], e$statistic)
  }
})

test_that("po_test tables, decision and notes", {
  r <- po_test(gs10 ~ tb3ms, term_structure)
  rs <- cointegrationtools:::po_rs[["2"]][["c"]][["Zt"]]
  expect_equal(r$critical_values,
               cointegrationtools:::rs_critical_values(rs, r$nobs + 1, 0))
  expect_equal(r$p.value, cointegrationtools:::rs_p_value(rs, r$statistic, r$nobs + 1, 0))
  ind <- simulate_coint(300, type = "independent", seed = 44)
  expect_false(po_test(y ~ x, ind)$reject_null)
  expect_match(po_test(y ~ x, ind, kernel = "qs")$warnings, "approximation", all = FALSE)
  expect_output(print(r), "Phillips-Ouliaris")
})

test_that("gh_test minimises over the break dates", {
  set.seed(45)
  n <- 120
  w <- cumsum(rnorm(n))
  d <- data.frame(x = w + rnorm(n, sd = 0.3))
  d$y <- ifelse(seq_len(n) <= 60, 1, 0.4) * d$x + rnorm(n, sd = 0.3)
  g <- gh_test(y ~ x, d, model = "C/S", lags = 1)
  i <- which.min(g$path$ADF)
  expect_equal(g$statistic, g$path$ADF[i])
  tb <- g$path$tb[i]
  D <- as.numeric(seq_len(n) > tb)
  u <- residuals(lm(y ~ D + x + I(x * D), d))
  expect_equal(g$statistic, unname(unitroottools::adf(u, trend = "n", lags = 1)$statistic))
  expect_equal(unname(residuals(g$regression)), unname(u))
  expect_equal(g$break_index, tb + 1)
  expect_true(g$reject_null)
  expect_true(abs(g$break_index - 61) <= 6)
  # Z statistics as po_stats on the same residuals
  j <- which.min(g$path$Za)
  Dj <- as.numeric(seq_len(n) > g$path$tb[j])
  uj <- residuals(lm(y ~ Dj + x + I(x * Dj), d))
  expect_equal(g$statistics$value[3], cointegrationtools:::po_stats(uj)$Za)
})

test_that("gh_test models, tables and checks", {
  d <- simulate_coint(150, beta = c(1, -0.5), seed = 46)
  for (md in c("C", "C/T", "C/S", "C/S/T")) {
    g <- gh_test(y ~ x1 + x2, d, model = md, lags = 0)
    tab <- cointegrationtools:::gh_tables[[md]]
    expect_equal(unname(g$critical_values), unname(tab$ADF[2, ]))
    expect_equal(g$statistics$`cv 5%`[3], tab$Za[2, 2])
    expect_equal(ncol(model.matrix(g$regression)),
                 switch(md, C = 4, `C/T` = 5, `C/S` = 6, `C/S/T` = 8))
  }
  expect_error(gh_test(y ~ x1 + x2, d, alpha = 0.025), "0.01, 0.05 or 0.10")
  expect_error(gh_test(y ~ x1 + x2, d, trim = 0.6), "trim")
  g <- gh_test(y ~ x1, d, model = "C", statistic = "Za", lags = 0)
  expect_equal(g$statistic, min(g$path$Za))
  expect_output(print(g), "asymptotic")
  expect_silent(plot(g))
  # the tables decrease with m and from 1% to 10% increase
  for (md in names(cointegrationtools:::gh_tables)) {
    tab <- cointegrationtools:::gh_tables[[md]]
    expect_true(all(diff(tab$ADF[, 2]) < 0) && all(apply(tab$ADF, 1, diff) > 0))
    expect_true(all(diff(tab$Za[, 2]) < 0) && all(apply(tab$Za, 1, diff) > 0))
  }
})

test_that("the fast ADF of gh_test() equals unitroottools::adf()", {
  for (s in 1:6) {
    set.seed(s)
    u <- as.numeric(arima.sim(list(ar = c(0.5, 0.3)), 150 + 30 * s))
    for (cr in c("BIC", "AIC")) {
      a <- unitroottools::adf(u, trend = "n", criterion = cr)
      f <- cointegrationtools:::adf_n_fast(u, NULL, a$lag_selection$kmax, cr)
      expect_equal(f$lags, a$lags)
      expect_equal(unname(f$statistic), unname(a$statistic))
    }
    f <- cointegrationtools:::adf_n_fast(u, 3L, NULL, "BIC")
    expect_equal(unname(f$statistic), unname(unitroottools::adf(u, trend = "n", lags = 3)$statistic))
  }
})
