\name{control.errorest}
\alias{control.errorest}
\title{ Control Error Rate Estimators }
\description{
  Some parameters that control the behaviour of \code{\link{errorest}}.
}
\usage{
control.errorest(k = 10, nboot = 25, strat = FALSE, random = TRUE, 
                 predictions = FALSE)
}
\arguments{
  \item{k}{integer, specify $k$ for $k$-fold cross-validation.}
  \item{nboot}{integer, number of bootstrap replications.}
  \item{strat}{logical, if \code{TRUE}, cross-validation is performed 
               using stratified sampling (for classification problems).}
  \item{random}{logical, if \code{TRUE}, cross-validation is performed using
                a random ordering of the data.}
  \item{predictions}{logical, indicates whether the prediction
                     for each observation should be returned or not
                     (classification and regression only). }
}
\value{
  A list with the same components as arguments. 
}
\keywords{misc}
