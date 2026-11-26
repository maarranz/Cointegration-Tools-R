sim_coint <- function(n = 250, seed = 2026) {
  set.seed(seed)
  w <- cumsum(rnorm(n))
  e <- rnorm(n, sd = 0.4)
  x1 <- w + e
  x2 <- cumsum(rnorm(n))
  # equilibrium error correlated with the innovations of x1 (endogeneity)
  z <- as.numeric(arima.sim(list(ar = 0.5), n, sd = 0.4)) + 0.5 * c(0, diff(w))
  data.frame(y = 2 + 1.5 * x1 - 0.7 * x2 + z, x1 = x1, x2 = x2)
}

det_matrix <- function(n, trend) {
  switch(trend, c = matrix(1, n, 1), ct = cbind(1, seq_len(n)))
}

test_that("IM-OLS reproduces cointReg exactly, with the scaled covariance", {
  skip_if_not_installed("cointReg")
  d <- sim_coint()
  X <- as.matrix(d[c("x1", "x2")])
  for (tr in c("c", "ct")) for (k in c("bartlett", "parzen", "qs")) {
    fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "IM", trend = tr, kernel = k))
    ref <- cointReg::cointRegIM(x = X, y = d$y, deter = det_matrix(nrow(d), tr),
                                kernel = substr(k, 1, 2))
    expect_equal(unname(coef(fit)), unname(ref$theta), tolerance = 1e-10)
    expect_equal(unname(sqrt(diag(fit$vcov))), unname(ref$sd.theta), tolerance = 1e-10)
  }
})

test_that("FM-OLS matches cointReg up to the T vs T - 1 bias factor", {
  skip_if_not_installed("cointReg")
  d <- sim_coint()
  X <- as.matrix(d[c("x1", "x2")])
  for (k in c("bartlett", "parzen", "qs")) {
    fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "FM", kernel = k))
    ref <- cointReg::cointRegFM(x = X, y = d$y, deter = det_matrix(nrow(d), "c"),
                                kernel = substr(k, 1, 2))
    expect_equal(fit$omega_uv, as.numeric(ref$omega.u.v), tolerance = 1e-10)
    expect_equal(unname(sqrt(diag(fit$vcov))), unname(ref$sd.theta), tolerance = 1e-10)
    expect_lt(max(abs(coef(fit) - ref$theta)), 2e-3)
  }
  fit_nw <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "FM", bandwidth = "nw"))
  ref_nw <- cointReg::cointRegFM(x = X, y = d$y, deter = det_matrix(nrow(d), "c"),
                                 bandwidth = "nw")
  expect_equal(fit_nw$omega_uv, as.numeric(ref_nw$omega.u.v), tolerance = 1e-10)
})

test_that("DOLS with given leads and lags reproduces cointReg and lm", {
  skip_if_not_installed("cointReg")
  d <- sim_coint()
  X <- as.matrix(d[c("x1", "x2")])
  fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "D", lags = 2, leads = 1))
  ref <- cointReg::cointRegD(x = X, y = d$y, deter = det_matrix(nrow(d), "c"),
                             n.lag = 2, n.lead = 1)
  expect_equal(unname(coef(fit)), unname(ref$theta), tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(fit$vcov))), unname(ref$sd.theta), tolerance = 1e-10)

  # the same regression written out with lm
  n <- nrow(d)
  dx1 <- c(NA, diff(d$x1)); dx2 <- c(NA, diff(d$x2))
  t <- 4:(n - 1)
  m <- lm(d$y[t] ~ d$x1[t] + d$x2[t] + dx1[t + 1] + dx2[t + 1] + dx1[t] + dx2[t] +
            dx1[t - 1] + dx2[t - 1] + dx1[t - 2] + dx2[t - 2])
  expect_equal(unname(coef(fit)), unname(coef(m)[1:3]), tolerance = 1e-10)
  expect_equal(nobs(fit), length(t))
})

test_that("DOLS chooses leads and lags on a common sample", {
  d <- sim_coint()
  fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "D", kmax = 3))
  sel <- fit$selection
  expect_equal(nrow(sel$table), 16)
  expect_equal(sel$nobs, nrow(d) - 1 - 2 * 3)
  expect_equal(c(fit$lags, fit$leads),
               unlist(sel$table[which.min(sel$table$value), c("lags", "leads")]),
               ignore_attr = TRUE)
  # the chosen pair refitted on its own sample equals a direct call
  direct <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "D",
                                       lags = fit$lags, leads = fit$leads))
  expect_equal(coef(fit), coef(direct))
})

test_that("residuals are the equilibrium error over the whole sample", {
  d <- sim_coint()
  for (m in c("OLS", "FM", "D", "IM")) {
    fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = m))
    b <- coef(fit)
    expect_equal(as.numeric(residuals(fit)),
                 d$y - b[["(Intercept)"]] - b[["x1"]] * d$x1 - b[["x2"]] * d$x2)
    expect_length(residuals(fit), nrow(d))
  }
  ols <- suppressMessages(coint_reg(y ~ x1 + x2, d))
  expect_equal(coef(ols), coef(lm(y ~ x1 + x2, d)))
})

test_that("input types give the same estimates", {
  d <- sim_coint()
  ref <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "FM"))
  tsd <- ts(as.matrix(d), start = c(1960, 1), frequency = 4)
  expect_equal(coef(coint_reg(y ~ x1 + x2, tsd, method = "FM")), coef(ref))
  expect_equal(coef(coint_reg(y ~ x1 + x2, zoo::as.zoo(tsd), method = "FM")), coef(ref))
  dd <- cbind(date = seq(as.Date("1960-01-01"), by = "quarter", length.out = nrow(d)), d)
  expect_equal(coef(coint_reg(y ~ x1 + x2, dd, method = "FM")), coef(ref))
  # sample window
  sub <- coint_reg(y ~ x1 + x2, tsd, start = c(1970, 1), end = c(2010, 4))
  expect_equal(sub$n, 164)
})

test_that("formula and argument errors are clear", {
  d <- sim_coint()
  expect_error(coint_reg(y ~ x1 - 1, d), "trend")
  expect_error(coint_reg(y ~ x1 * x2, d), "Interactions")
  expect_error(suppressMessages(coint_reg(log(y - 100) ~ x1, d)), "missing or infinite")
  expect_error(coint_reg(y ~ x1, d, method = "D", lags = -1), "non-negative")
  expect_error(suppressMessages(coint_reg(y ~ z, d)), "not found")
})

test_that("OLS covariance warns; Wald tests need FM, D or IM", {
  d <- sim_coint()
  ols <- suppressMessages(coint_reg(y ~ x1 + x2, d))
  expect_warning(vcov(ols), "not valid")
  expect_error(coint_wald_test(ols, "x1 = 1.5"), "FM")
})

test_that("Wald test: text and matrix forms agree with the formula", {
  d <- sim_coint()
  fit <- suppressMessages(coint_reg(y ~ x1 + x2, d, method = "D"))
  b <- coef(fit); V <- fit$vcov
  w1 <- coint_wald_test(fit, c("x1 = 1.5", "x1 + x2 = 0.8"))
  R <- rbind(c(0, 1, 0), c(0, 1, 1)); r <- c(1.5, 0.8)
  dd <- R %*% b - r
  expect_equal(w1$statistic, as.numeric(t(dd) %*% solve(R %*% V %*% t(R), dd)))
  expect_equal(w1$df, 2)
  w2 <- coint_wald_test(fit, R, r)
  expect_equal(w2$statistic, w1$statistic)
  expect_equal(coint_wald_test(fit, "x1 - 1.5 = 0")$statistic,
               coint_wald_test(fit, "x1 = 1.5")$statistic)
  expect_equal(coint_wald_test(fit, "`(Intercept)` = 2")$df, 1)
  expect_error(coint_wald_test(fit, "x1 * x2 = 1"), "not a linear")
  expect_error(coint_wald_test(fit, "x3 = 1"), "unknown coefficient")
  expect_output(print(w1), "Chi-square\\(2\\)")
})
