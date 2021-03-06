% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{sigma.blblm}
\alias{sigma.blblm}
\title{Calculate the residual standard error estimate from the Bootstrapped estimates
Construct a confidence interval for the sigma estimate if specified}
\usage{
\method{sigma}{blblm}(object, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{blblm fitted object}

\item{confidence}{logical Indicate whether you would like a confidence interval for sigma.}

\item{level}{numeric Confidence level.}

\item{...}{further arguments}
}
\value{
numeric Sigma estimate or sigma estimate with lower and upper bounds of confidence interval.
}
\description{
Calculate the residual standard error estimate from the Bootstrapped estimates
Construct a confidence interval for the sigma estimate if specified
}
\examples{
sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), confidence = TRUE, level = 0.95)
}
