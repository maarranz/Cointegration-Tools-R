test_that("integration_order classifies I(0), I(1) and I(2) series", {
  set.seed(42)
  n <- 300
  d <- data.frame(u = rnorm(n), w = cumsum(rnorm(n)))
  d$w2 <- cumsum(cumsum(rnorm(n)))
  io <- integration_order(~ u + w + w2, d)
  expect_s3_class(io, "integration_order")
  expect_equal(io$variable, c("u", "w", "w2"))
  expect_equal(io$order, c("I(0)", "I(1)", "I(2) or higher"))
  expect_named(attr(io, "strategies"), c("u", "w", "w2"))
  expect_output(print(io), "Order of integration")
  # a two-sided formula and transformations work too
  d$p <- exp(d$w / 10)
  io2 <- integration_order(log(p) ~ w, d)
  expect_equal(io2$variable, c("log(p)", "w"))
  expect_equal(io2$level_stat[1], io$level_stat[2], tolerance = 1e-8)
})

test_that("integration_order: hamilton strategy and argument checks", {
  set.seed(1)
  d <- data.frame(a = cumsum(rnorm(200)), b = cumsum(rnorm(200)))
  io <- integration_order(~ a + b, d, strategy = "hamilton", trending = c(FALSE, FALSE))
  expect_equal(io$order, c("I(1)", "I(1)"))
  expect_error(integration_order(~ a + b, d, trending = c(TRUE, FALSE, TRUE)), "one per variable")
  expect_error(integration_order("a", d), "formula")
})

test_that("simulate_coint returns the stated DGPs", {
  d <- simulate_coint(300, beta = c(1.5, -0.7), alpha = 2, rho = 0.3, seed = 7)
  expect_named(d, c("t", "y", "x1", "x2", "z"))
  expect_equal(d$y, 2 + 1.5 * d$x1 - 0.7 * d$x2 + d$z)
  expect_equal(d$z[-1] - 0.3 * d$z[-300], d$z[-1] - 0.3 * d$z[-300])  # AR(1) by construction
  expect_equal(simulate_coint(50, seed = 3), simulate_coint(50, seed = 3))
  ind <- simulate_coint(100, type = "independent", seed = 1)
  expect_false("z" %in% names(ind))
  br <- simulate_coint(200, beta = 1, beta_post = 0.4, type = "break", seed = 2)
  expect_equal(br$y[150], 0.4 * br$x[150] + br$z[150])
  expect_equal(br$y[50], br$x[50] + br$z[50])
  # endogeneity: correlation between dx and the error innovation
  en <- simulate_coint(5000, rho = 0, corr = 0.6, seed = 4)
  expect_equal(cor(diff(en$x), en$z[-1]), 0.6, tolerance = 0.05)
  expect_error(simulate_coint(100, rho = 1), "stationary")
})

test_that("simulate_spurious: rejections grow with T", {
  s <- simulate_spurious(T = c(50, 400), M = 300, seed = 1)
  expect_s3_class(s, "simulate_spurious")
  expect_true(all(s$reject > 0.3))
  expect_gt(s$reject[2], s$reject[1])
  expect_equal(nrow(attr(s, "draws")), 600)
  # the t ratio equals lm's
  set.seed(5); y <- cumsum(rnorm(60)); x <- cumsum(rnorm(60))
  set.seed(5); s1 <- simulate_spurious(T = 60, M = 1)
  expect_equal(attr(s1, "draws")$t, coef(summary(lm(y ~ x)))["x", "t value"])
  expect_output(print(s), "share")
})

test_that("term_structure is documented data", {
  expect_true(all(c("date", "tb3ms", "gs10") %in% names(term_structure)))
  expect_s3_class(term_structure$date, "Date")
})
