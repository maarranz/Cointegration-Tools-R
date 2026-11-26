sim_pair <- function(n = 200, seed = 1, coint = TRUE) {
  set.seed(seed)
  w <- cumsum(rnorm(n))
  d <- data.frame(x = w + rnorm(n, sd = 0.4), x2 = cumsum(rnorm(n)))
  d$y <- if (coint) 2 + 1.5 * d$x + as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4))
         else cumsum(rnorm(n))
  d
}

test_that("the statistic is the ADF t ratio of the Step 1 residuals, no deterministics", {
  d <- sim_pair()
  for (tr in c("c", "ct", "n")) {
    r <- suppressMessages(eg_test(y ~ x + x2, d, trend = tr, lags = 2))
    D <- switch(tr, c = "1", ct = "1 + seq_along(y)", n = "0")
    z <- residuals(lm(as.formula(paste("y ~", D, "+ x + x2")), d))
    a <- unitroottools::adf(z, trend = "n", lags = 2)
    expect_equal(r$statistic, a$statistic)
    expect_equal(r$nobs, a$nobs)
    expect_equal(r$nvars, 3L)
  }
})

test_that("critical values and p-values come from the Engle-Granger surfaces", {
  d <- sim_pair()
  r <- suppressMessages(eg_test(y ~ x, d, lags = 0))
  rs <- cointegrationtools:::eg_rs[["2"]][["c"]]
  T_eq <- r$nobs + 1L
  expect_equal(r$critical_values, cointegrationtools:::rs_critical_values(rs, T_eq, 0))
  # p-value equals the level at the corresponding critical value
  for (a in c(0.01, 0.05, 0.10)) {
    expect_equal(cointegrationtools:::rs_p_value(rs, r$critical_values[[as.character(a)]], T_eq, 0),
                 a, tolerance = 1e-8)
  }
  # more negative than the plain ADF critical values with a constant
  adf_cv <- unitroottools::adf(d$x, trend = "c", lags = 0)$critical_values
  expect_true(all(r$critical_values < adf_cv))
})

test_that("asymptotic values are close to MacKinnon (2010)", {
  # MacKinnon (2010), Table 2, N = 2, constant: tau_inf at 1/5/10%
  rs <- cointegrationtools:::eg_rs[["2"]][["c"]]
  q <- cointegrationtools:::rs_critical_values(rs, 1e8, 0)
  expect_equal(unname(q), c(-3.89644, -3.33613, -3.04445), tolerance = 0.01)
})

test_that("decisions: cointegrated pair rejects, independent walks do not", {
  expect_true(suppressMessages(eg_test(y ~ x, sim_pair()))$reject_null)
  # with this sample MAIC picks 14 lags and loses the rejection (why BIC is the default)
  expect_false(suppressMessages(eg_test(y ~ x, sim_pair(), criterion = "MAIC"))$reject_null)
  expect_false(suppressMessages(eg_test(y ~ x, sim_pair(coint = FALSE)))$reject_null)
})

test_that("normalisations repeat the test with each variable on the left", {
  d <- sim_pair(coint = FALSE)
  r <- suppressMessages(eg_test(y ~ x + x2, d, normalisations = TRUE))
  expect_equal(r$normalisations$dependent, c("y", "x", "x2"))
  expect_equal(r$normalisations$statistic[1], r$statistic)
  r2 <- suppressMessages(eg_test(x ~ y + x2, d))
  expect_equal(r$normalisations$statistic[2], r2$statistic)
  expect_output(print(r), "Each variable as the dependent one")
})

test_that("too many variables and printing", {
  set.seed(3)
  d <- as.data.frame(matrix(cumsum(rnorm(700)), 100, 7))
  expect_error(suppressMessages(eg_test(V1 ~ V2 + V3 + V4 + V5 + V6 + V7, d)), "2 to 6")
  d <- sim_pair()
  d$y <- 1 + d$x + rnorm(nrow(d), sd = 0.1)        # white-noise equilibrium error
  r <- suppressMessages(eg_test(y ~ x, d))
  expect_output(print(r), "not the ADF tables")
  expect_output(print(r), "< 0.001", fixed = TRUE)
  expect_true(any(grepl("beyond the simulated quantiles", r$warnings)))
})
