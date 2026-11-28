test_that("the ECM statistic is the plain-OLS t ratio of y(-1)", {
  d <- simulate_coint(200, beta = c(1, -0.5), rho = 0.7, seed = 31)
  for (tr in c("c", "ct")) {
    r <- suppressMessages(ecm_test(y ~ x1 + x2, d, trend = tr, p = 2))
    n <- 200; rows <- 4:n
    dy <- c(NA, diff(d$y)); d1 <- c(NA, diff(d$x1)); d2 <- c(NA, diff(d$x2))
    df <- data.frame(dy = dy[rows], yl = d$y[rows - 1], x1l = d$x1[rows - 1], x2l = d$x2[rows - 1],
                     d1 = d1[rows], d2 = d2[rows], dyl1 = dy[rows - 1], d1l1 = d1[rows - 1],
                     d2l1 = d2[rows - 1], dyl2 = dy[rows - 2], d1l2 = d1[rows - 2],
                     d2l2 = d2[rows - 2], tt = rows)
    m <- if (tr == "c") lm(dy ~ . - tt, df) else lm(dy ~ ., df)
    expect_equal(r$statistic, coef(summary(m))["yl", "t value"])
    expect_equal(unname(r$long_run), unname(-coef(m)[c("x1l", "x2l")] / coef(m)[["yl"]]))
    expect_equal(r$nobs, length(rows))
  }
})

test_that("the 'unit' parameterisation of the chapter gives the same t ratio", {
  d <- simulate_coint(150, rho = 0.5, seed = 32)
  r <- suppressMessages(ecm_test(y ~ x, d, p = 0))
  dy <- diff(d$y); dx <- diff(d$x); yl <- d$y[-150]; xl <- d$x[-150]
  m <- lm(dy ~ dx + I(yl - xl) + xl)
  expect_equal(r$statistic, coef(summary(m))["I(yl - xl)", "t value"])
})

test_that("tables, decisions and lag selection", {
  r <- suppressMessages(ecm_test(gs10 ~ tb3ms, term_structure))
  rs <- cointegrationtools:::ecm_rs[["2"]][["c"]]
  expect_equal(r$critical_values,
               cointegrationtools:::rs_critical_values(rs, r$nobs + r$p + 1, r$p))
  expect_true(r$reject_null)
  expect_equal(nrow(r$selection), 5)
  expect_equal(r$p, r$selection$p[which.min(r$selection$value)])
  ind <- simulate_coint(300, type = "independent", seed = 33)
  expect_false(suppressMessages(ecm_test(y ~ x, ind))$reject_null)
  # the ECM 5% value is less negative than Engle-Granger's for n = 2
  eg <- suppressMessages(eg_test(y ~ x, ind, lags = 0))
  expect_gt(r$critical_values[["0.05"]], eg$critical_values[["0.05"]])
})

test_that("argument checks and printing", {
  d <- simulate_coint(100, seed = 34)
  expect_error(suppressMessages(ecm_test(y ~ x, d, p = 5)), "p = 0 to 4")
  expect_error(suppressMessages(ecm_test(y ~ x, d, trend = "n")), "should be one of")
  r <- suppressMessages(ecm_test(y ~ x, d, p = 1))
  expect_output(print(r), "weakly exogenous")
  expect_output(print(r), "Long run")
})
