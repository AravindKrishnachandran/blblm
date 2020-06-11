test_that("blblm works", {

  fit <- blblm(mpg ~ wt * hp, data = mtcars, parallel = FALSE, m = 3, B = 100)
  expect_s3_class(fit, "blblm")

  expect_equal(dim(predict(fit, data.frame(wt = c(2.5, 3),
                                           hp = c(150, 170)), confidence = TRUE)),
        c(2,3))


  expect_length(confint(fit, c("wt", "hp", "wt:hp")), 6)

  expect_equal(parallelization(TRUE, workers = 2), TRUE)

  expect_length(sigma(fit, confidence = TRUE), 3)

})
