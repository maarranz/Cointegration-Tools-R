test_that("the package loads and has its help page", {
  expect_true("cointegrationtools" %in% loadedNamespaces())
  expect_true(nzchar(utils::packageDescription("cointegrationtools")$Title))
})
