#' @import purrr
#' @import stats
#' @import furrr
#' @import vroom
#' @import future
#' @import nycflights13
#' @importFrom utils capture.output
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Gives the user the choice if he/she would like to use parallelization in the blblm calculations.
#' User can specify number of workers.
#'
#' @param parallel logical Indicate whether or not you want to use parallelization.
#' @param workers numeric Specify number of workers to distribute tasks.
#'
#' @return logical TRUE if using parallelization, FALSE if not.
#'
#' @export
#'
#' @examples
#' parallelization(FALSE)
parallelization <- function(parallel = FALSE, workers) {
  if (parallel) {
    suppressWarnings(plan(multiprocess, workers = workers))
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Compute Linear Regression Model using Bag of Little Bootstraps algorithm.
#' The user can specify whether or not she would like to use parallelization.
#' Parallelization is done using future_map.
#' If a character vector of files is passed in, each file is read in the workers.
#' If an entire data frame is passed in, the data is split and it will either use future_map or regular map to calculate coefficents and sigma.
#'
#' @param formula Specify formula for model.
#' @param data character or data.frame Character vector listing files or whole dataset.
#' @param parallel logical TRUE if you specified TRUE in parallelization function, FALSE otherwise.
#' @param m numeric Number of subsamples to split the data if whole dataset is specified.
#' @param B numeric Number of bootstrap estimates.
#'
#' @return class blblm object Can extract coefficent, sigma estimates, along with formula.
#'
#' @export
#'
#' @examples
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
blblm <- function(formula, data, parallel = FALSE, m = 10, B = 5000) {
  if (parallel) {
    if (is.character(data)) {
      n <- length(vroom_lines(data, altrep = TRUE, progress = FALSE)) - length(data)
      estimates <- data %>% future_map(~ {
        lm_each_subsample(formula = formula, data = vroom(., col_types = cols()), n = n, B = B)
      })
      res <- list(estimates = estimates, formula = formula)
      class(res) <- "blblm"
      invisible(res)
    } else {
      estimates <- future_map(
        split_data(data, m),
        ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
      )
      res <- list(estimates = estimates, formula = formula)
      class(res) <- "blblm"
      invisible(res)
    }
  } else {
    data_list <- split_data(data, m)
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
    )
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}


#' Split data into m parts of approximated equal sizes
#'
#' @param data data.frame
#' @param m numeric
#'
#' @return subsamples of data
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' Compute the estimates by repeating lm_each_boot B times for a given subsample.
#'
#' @param formula formula
#' @param data data.frame
#' @param n numeric
#' @param B numeric
#'
#' @return list of coefficents and sigma
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' Compute the regression estimates for a blb dataset by including frequencies.
#'
#' @param formula formula
#' @param data data.frame
#' @param n numeric
#'
#' @return list of coefficents and sigma
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' Estimate the regression estimates based on given number of repetitions
#'
#' @param formula formula
#' @param data data.frame
#' @param freqs matrix
#'
#' @return list of coefficents and sigma
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' Compute the coefficients from fit
#'
#' @param fit lm, glm
#'
#' @return numeric coefficent estimates
blbcoef <- function(fit) {
  coef(fit)
}


#' Compute sigma (residual standard error) from fit
#'
#' @param fit lm
#'
#' @return numeric sigma estimate
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the fitted blblm formula
#'
#' @param x blblm fitted model
#' @param ... further arguments
#'
#' @return character The Model Formula.
#'
#' @export
#'
#' @method print blblm
#'
#' @examples
#' print(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Calculate the residual standard error estimate from the Bootstrapped estimates
#' Construct a confidence interval for the sigma estimate if specified
#'
#' @param object blblm fitted object
#' @param confidence logical Indicate whether you would like a confidence interval for sigma.
#' @param level numeric Confidence level.
#' @param ... further arguments
#'
#' @return numeric Sigma estimate or sigma estimate with lower and upper bounds of confidence interval.
#'
#' @export
#'
#' @method sigma blblm
#'
#' @examples
#' sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), confidence = TRUE, level = 0.95)
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Get the estimated values for the coefficients from the model.
#'
#' @param object blblm fitted object
#' @param ... further arguments
#'
#' @return numeric The average of the Bootstrapped coefficients from each subsample.
#'
#' @export
#'
#' @method coef blblm
#'
#' @examples
#' coef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Construct a confidence interval for named coefficients in the model.
#' If coefficient names not specified, than confidence intervals constructed for all coeffienct names.
#'
#' @param object blblm fitted model
#' @param parm character Names of the coeffiecients used for confidence intervals.
#' @param level numeric Confidence level.
#' @param ... further arguments
#'
#' @return matrix Confidence Interval lower and upper bounds for each coefficient.
#'
#' @export
#'
#' @method confint blblm
#'
#' @examples
#' confint(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), c("wt", "hp"))
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Predict a response for new data points using model formula.
#' Also calculate confidence intervale for these response estimates if specified.
#'
#' @param object blblm fitted model
#' @param new_data numeric New data points to predict.
#' @param confidence logical Specify whether or not you would like a confidence interval.
#' @param level numeric Confidence level.
#' @param ... further arguments
#'
#' @return numeric Predicted estimates with or without a confidence interval.
#'
#' @export
#'
#' @method predict blblm
#'
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
#' predict(fit, data.frame(wt = c(7, 5.5), hp = c(130, 190)), confidence = TRUE, level = 0.95)
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms.formula(object$formula, data = new_data), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

#' Calculates the mean of an object and the lower and upper bounds for it
#'
#' @param x numeric Vector of estimates.
#' @param level numeric Confidence level.
#'
#' @return numeric estimate with bounds
mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

#' Applies a function to an object via map and then calculcates the mean of the result
#'
#' @param .x object
#' @param .f function
#' @param ... further arguments
#'
#' @return numeric mean
map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

#' Applies a function to an object via map and binds the columns
#'
#' @param .x object
#' @param .f function
#' @param ... further arguments
#'
#' @return tibble
map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

#' Applies a function to an object via map and binds the rows
#'
#' @param .x object
#' @param .f function
#' @param ... further arguments
#'
#' @return tibble
map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
