\name{control.errorest}
\alias{control.errorest}
\title{ Control Error Rate Estimators }
\description{
  Some parameters that control the behaviour of \code{\link{errorest}}.
}
\usage{
control.errorest(k = 10, nboot = 25, strat = FALSE, random = TRUE, 
                 predictions = FALSE, getmodels=FALSE, list.tindx = NULL)
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
                     (classification and regression only). For a bootstrap
based estimator a matrix of size 'number of observations' times nboot is
returned with predicted values of the ith out-of-bootstrap sample in column
i and 'NA's for those observations not included in the ith out-of-bootstrap 
sample.}
  \item{getmodels}{logical, indicates a list of all models should be
    returned. For cross-validation only.}
  \item{list.tindx}{list of numeric vectors, indicating which
    observations are included in each bootstrap or cross-validation sample, respectively.}
}
\value{
  A list with the same components as arguments. 
}
\keyword{misc}
