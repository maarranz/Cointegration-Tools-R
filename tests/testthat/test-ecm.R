test_that("coint_ecm equals the two-step regression written out with lm", {
  d <- simulate_coint(250, beta = 1.5, rho = 0.6, seed = 21)
  lr <- suppressMessages(coint_reg(y ~ x, d, method = "FM"))
  ecm <- coint_ecm(lr, p = 2, q = 1)
  expect_s3_class(ecm, c("coint_ecm", "adl"))
  z <- as.numeric(residuals(lr))
  dy <- c(NA, diff(d$y)); dx <- c(NA, diff(d$x)); ec <- c(NA, z[-250])
  L <- function(v, j) c(rep(NA, j), v[1:(250 - j)])
  m <- lm(dy ~ L(dy, 1) + L(dy, 2) + dx + L(dx, 1) + ec)
  expect_equal(unname(coef(ecm)[sort(names(coef(ecm)))]),
               unname(coef(m)[c("(Intercept)", "dx", "L(dx, 1)", "L(dy, 1)", "L(dy, 2)", "ec")]),
               tolerance = 1e-10)
  expect_equal(ecm$adjustment[["estimate"]], coef(m)[["ec"]])
  expect_lt(ecm$adjustment[["estimate"]], 0)
  expect_output(print(ecm), "Adjustment")
})

test_that("adltools diagnostics and stability tests run on coint_ecm", {
  lr <- suppressMessages(coint_reg(gs10 ~ tb3ms, term_structure, method = "FM"))
  ecm <- coint_ecm(lr, p = 1, q = 1)
  expect_s3_class(adltools::cusum_test(ecm), "adltools_test")
  expect_no_error(adltools::adl_diagnostics(ecm))
  expect_lt(ecm$adjustment[["estimate"]], 0)
})

test_that("weak exogeneity test equals its marginal regression", {
  d <- simulate_coint(250, beta = 1, rho = 0.5, seed = 22)
  lr <- suppressMessages(coint_reg(y ~ x, d, method = "D"))
  w <- coint_weak_exogeneity_test(lr, lags = 2)
  z <- as.numeric(residuals(lr))
  dy <- c(NA, diff(d$y)); dx <- c(NA, diff(d$x)); ec <- c(NA, z[-250])
  L <- function(v, j) c(rep(NA, j), v[1:(250 - j)])
  m <- lm(dx ~ ec + L(dy, 1) + L(dy, 2) + L(dx, 1) + L(dx, 2))
  expect_equal(w$kappa, coef(m)[["ec"]])
  expect_equal(w$statistic, coef(summary(m))["ec", "t value"])
  expect_identical(coint_weak_exogeneity_test(coint_ecm(lr))$kappa,
                   coint_weak_exogeneity_test(lr)$kappa)
  # x is a random walk unaffected by z in simulate_coint: rarely rejected
  expect_gt(w$p.value, 0.01)
})

test_that("simulate_eg reproduces the EG tables and separates the three tests", {
  s <- simulate_eg(T = 200, M = 2000, seed = 3)
  expect_s3_class(s, "simulate_eg")
  s0 <- as.data.frame(unclass(s), check.names = FALSE)
  eg <- unlist(s0[s0$statistic == "EG", -1]); tabs <- unlist(s0[s0$statistic == "EG tables", -1])
  expect_equal(unname(eg), unname(tabs), tolerance = 0.12)
  adf <- unlist(s0[s0$statistic == "ADF", -1])
  expect_true(all(eg < adf))
  expect_output(print(s), "Null distributions")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(s))
  expect_error(simulate_eg(n = 9), "between 2 and 6")
})
