test_that("the last recursive estimate is the whole-sample fit", {
  d <- simulate_coint(200, beta = 1.5, rho = 0.5, seed = 11)
  for (m in c("FM", "IM", "OLS", "D")) {
    rec <- suppressMessages(coint_recursive(y ~ x, d, method = m))
    full <- attr(rec, "full")
    last <- rec[rec$end == max(rec$end), ]
    expect_equal(setNames(last$estimate, last$term), coef(full))
    expect_equal(min(rec$end), max(30, ceiling(0.15 * 200)))
  }
})

test_that("rolling windows equal coint_reg on the same window", {
  d <- simulate_coint(160, beta = 1, rho = 0.3, seed = 12)
  d$t <- NULL
  roll <- suppressMessages(coint_recursive(y ~ x, d, method = "FM", window = 80, trend = "ct"))
  one <- roll[roll$end == 120, ]
  ref <- suppressMessages(coint_reg(y ~ x, d, method = "FM", trend = "ct", start = 41, end = 120))
  expect_equal(setNames(one$estimate, one$term), coef(ref))
  expect_equal(one$std.error, unname(sqrt(diag(ref$vcov))))
  expect_equal(length(unique(roll$end)), 160 - 80 + 1)
})

test_that("DOLS keeps the whole-sample leads and lags; arguments are checked", {
  d <- simulate_coint(200, seed = 13)
  rec <- suppressMessages(coint_recursive(y ~ x, d, method = "D", kmax = 2))
  full <- attr(rec, "full")
  r100 <- rec[rec$end == 100, ]
  ref <- suppressMessages(coint_reg(y ~ x, d, method = "D", lags = full$lags, leads = full$leads,
                                    end = 100))
  expect_equal(setNames(r100$estimate, r100$term), coef(ref))
  expect_error(suppressMessages(coint_recursive(y ~ x, d, window = 5)), "between")
  expect_error(suppressMessages(coint_recursive(y ~ x, d, min_obs = 500)), "smaller")
})

test_that("print and plot work", {
  d <- simulate_coint(120, seed = 14)
  rec <- suppressMessages(coint_recursive(y ~ x, d))
  expect_output(print(rec), "Recursive estimates")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(rec, ref = c(x = 1)))
  expect_error(plot(rec, terms = "z"), "Unknown")
})
