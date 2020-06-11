test_that("blblm works", {

  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_s3_class(fit, "blblm")


  expect_length(confint(fit, c("wt", "hp", "wt:hp")), 6)

  expect_equal(parallelization(TRUE, workers = 2), TRUE)

  expect_length(sigma(fit, confidence = TRUE), 3)

})
