# Regression tests for the findings of the audit of 26 Nov 2026 (H = high,
# M = medium, L = low), built from the audit's own counterexamples.

test_that("H01: an exact Step 1 fit is refused, whatever the scale", {
  set.seed(1)
  d <- data.frame(x = cumsum(rnorm(100)))
  d$y <- 2 * d$x
  expect_error(suppressMessages(eg_test(y ~ x, d, lags = 0)), "fits exactly")
  d$y <- 1e6 + 2e-3 * d$x
  expect_error(suppressMessages(eg_test(y ~ x, d, lags = 0)), "fits exactly")
  d$y <- 2 * d$x
  expect_error(suppressMessages(eg_test(y ~ x, d, trend = "n", lags = 0)), "fits exactly")
})

test_that("H02: with one order fixed, DOLS selects the other given it", {
  set.seed(33)
  d <- simulate_coint(150, beta = 1, rho = 0.5, corr = 0.5)
  fit <- suppressMessages(coint_reg(y ~ x, d, method = "D", lags = 3, kmax = 3))
  expect_equal(fit$lags, 3)
  tab <- fit$selection$table
  expect_true(all(tab$lags == 3))
  expect_equal(fit$leads, tab$leads[which.min(tab$value)])
  # independent grid: leads 0..3 with lags 3, on rows 2+3 .. T-3
  n <- nrow(d); rows <- (2 + 3):(n - 3)
  dx <- c(NA, diff(d$x))
  ic <- sapply(0:3, function(q) {
    Z <- cbind(1, d$x[rows], sapply(-q:3, function(j) dx[rows - j]))
    e <- lm.fit(Z, d$y[rows])$residuals
    length(rows) * log(sum(e^2) / length(rows)) + 2 * ncol(Z)
  })
  expect_equal(fit$leads, which.min(ic) - 1)
  fit2 <- suppressMessages(coint_reg(y ~ x, d, method = "D", leads = 1, kmax = 3))
  expect_true(all(fit2$selection$table$leads == 1))
  expect_output(print(fit2), "leads fixed")
})

test_that("H03: only linear restrictions are accepted, without evaluating code", {
  d <- simulate_coint(200, beta = 1.5, seed = 2)
  fit <- suppressMessages(coint_reg(y ~ x, d, method = "FM"))
  expect_error(coint_wald_test(fit, "x + x*(x-1)*(x-1.79) = 1"), "product of coefficients")
  expect_error(coint_wald_test(fit, "exp(x) = 1"), "not allowed")
  expect_error(coint_wald_test(fit, "x = system('echo hi')"), "not allowed")
  expect_error(coint_wald_test(fit, "x / x = 1"), "division by a coefficient")
  expect_equal(coint_wald_test(fit, "2*x/2 - 1.5 = 0")$statistic,
               coint_wald_test(fit, "x = 1.5")$statistic)
  expect_equal(coint_wald_test(fit, "-(1.5 - x) = 0")$statistic,
               coint_wald_test(fit, "x = 1.5")$statistic)
})

test_that("H04: text dates must be ISO dates, and gaps are then detected", {
  d <- simulate_coint(120, seed = 3)
  d$date <- format(seq(as.Date("2000-01-01"), by = "month", length.out = 120))
  expect_equal(coef(coint_reg(y ~ x, d, method = "FM", time = "date")),
               coef(suppressMessages(coint_reg(y ~ x, d[names(d) != "date"], method = "FM"))))
  expect_error(coint_reg(y ~ x, d[-35, ], method = "FM", time = "date"), "gaps")
  d$date <- sprintf("t%03d", 1:120)
  expect_error(coint_reg(y ~ x, d, time = "date"), "ISO dates")
})

test_that("H05: the EG result states the I(1) assumption", {
  set.seed(22)
  d <- data.frame(x = cumsum(rnorm(200)), y = rnorm(200))
  r <- suppressMessages(eg_test(y ~ x, d))
  expect_match(r$assumption, "I\\(1\\)")
  expect_output(print(r), "Assumes:")
})

test_that("M01: infeasible DOLS orders and non-finite counts are refused", {
  d <- simulate_coint(30, seed = 4)
  expect_error(suppressMessages(coint_reg(y ~ x, d, method = "D", lags = 20, leads = 0)),
               "needs more observations")
  expect_error(suppressMessages(coint_reg(y ~ x, d, method = "D", kmax = 12)),
               "needs more observations")
  expect_error(coint_reg(y ~ x, d, method = "D", lags = Inf), "non-negative")
})

test_that("M02: levels are validated and named R columns matched by name", {
  d <- simulate_coint(200, beta = 1.5, seed = 5)
  fit <- suppressMessages(coint_reg(y ~ x, d, method = "FM"))
  expect_error(coint_wald_test(fit, "x = 1", alpha = 2), "between 0 and 1")
  expect_error(summary(fit, level = 95), "between 0 and 1")
  expect_error(suppressMessages(eg_test(y ~ x, d, alpha = 0)), "between 0 and 1")
  R1 <- matrix(c(0, 1), 1, dimnames = list(NULL, c("(Intercept)", "x")))
  R2 <- matrix(c(1, 0), 1, dimnames = list(NULL, c("x", "(Intercept)")))
  expect_equal(coint_wald_test(fit, R1, 1.5)$statistic, coint_wald_test(fit, R2, 1.5)$statistic)
  expect_equal(coint_wald_test(fit, R1, 1.5)$statistic, coint_wald_test(fit, "x = 1.5")$statistic)
  R3 <- matrix(c(1, 0), 1, dimnames = list(NULL, c("a", "b")))
  expect_error(coint_wald_test(fit, R3, 1), "column names")
  expect_error(coint_wald_test(fit, matrix(c(0, NA), 1), 1), "finite")
})

test_that("M03: p equals the level at the tabulated critical values", {
  rs <- cointegrationtools:::eg_rs[["3"]][["ct"]]
  for (T in c(40, 100, 333)) for (k in c(0, 4)) {
    cv <- cointegrationtools:::rs_critical_values(rs, T, k)
    for (a in c(0.01, 0.05, 0.10)) {
      p <- cointegrationtools:::rs_p_value(rs, cv[[format(a)]], T, k)
      expect_identical(p, rs$probs[which.min(abs(rs$probs - a))])
    }
  }
  # a bounded p with alpha below the bound gives no decision
  d <- simulate_coint(200, beta = 1, rho = 0, sd_z = 0.1, seed = 6)
  r <- suppressMessages(eg_test(y ~ x, d, alpha = 0.0005))
  expect_true(is.na(r$reject_null))
  expect_output(print(r), "undetermined")
})

test_that("M04: local objects in formulas; regressor names that clash", {
  d <- simulate_coint(150, seed = 7)
  sh <- function(v) v + 100
  k <- 2
  f1 <- suppressMessages(coint_reg(y ~ sh(x), d, method = "FM"))
  f2 <- suppressMessages(coint_reg(y ~ I(k * x), d, method = "FM"))
  expect_equal(unname(coef(f1)[2]), unname(coef(suppressMessages(coint_reg(y ~ x, d, method = "FM")))[2]))
  # rescaling x changes the Andrews bandwidth slightly (equal column weights)
  expect_equal(unname(coef(f2)[2]) * 2, unname(coef(f1)[2]), tolerance = 1e-2)
  d$trend <- d$x
  expect_error(coint_reg(y ~ trend, d, trend = "ct"), "rename")
})

test_that("M05: simulator arguments are validated", {
  expect_error(simulate_coint(100, sd_x = 0), "positive")
  expect_error(simulate_coint(100, beta = NA), "finite")
  expect_error(simulate_spurious(T = 2, M = 10), "at least 5")
  expect_error(simulate_spurious(T = 50, M = 0), "at least 1")
})

test_that("L01: ordinary p-values are printed as computed", {
  d <- simulate_coint(300, beta = 1.5, seed = 8)
  fit <- suppressMessages(coint_reg(y ~ x, d, method = "FM"))
  w <- coint_wald_test(fit, "x = 0")
  expect_lt(w$p.value, 1e-10)
  expect_false(any(grepl("< 0.001", capture.output(print(w)), fixed = TRUE)))
})
