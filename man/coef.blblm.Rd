% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{coef.blblm}
\alias{coef.blblm}
\title{Get the estimated values for the coefficients from the model.}
\usage{
\method{coef}{blblm}(object, ...)
}
\arguments{
\item{object}{blblm fitted object}

\item{...}{further arguments}
}
\value{
numeric The average of the Bootstrapped coefficients from each subsample.
}
\description{
Get the estimated values for the coefficients from the model.
}
\examples{
coef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
}
