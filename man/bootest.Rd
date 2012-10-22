\name{bootest}
\alias{bootest}
\alias{bootest.default}
\alias{bootest.factor}
\alias{bootest.numeric}
\alias{bootest.integer}
\alias{bootest.Surv}
\title{Bootstrap Error Rate Estimators}
\description{
  Those functions are low-level functions used by \code{\link{errorest}} and
are normally not called by users.
}
\usage{
\method{bootest}{factor}(y, formula, data, model, predict, nboot=25,
bc632plus=FALSE, list.tindx = NULL, predictions = FALSE, 
both.boot = FALSE, \dots)}
\arguments{
  \item{y}{the response variable, either of class \code{factor}
(classification), \code{numeric} (regression) or \code{Surv} (survival).}
  \item{formula}{a formula object.}
  \item{data}{data frame of predictors and response described in
\code{formula}.}   
  \item{model}{a function implementing the predictive model to be
evaluated. The function \code{model} can either return an
               object representing a fitted model or a function with
               argument \code{newdata} which returns predicted values. In
               this case, the \code{predict} argument to \code{errorest} is
               ignored.}
  \item{predict}{a function with arguments \code{object} and \code{newdata}
only which predicts the status of the observations in \code{newdata} based
on the fitted model in \code{object}.}
  \item{nboot}{number of bootstrap replications to be used.}
  \item{bc632plus}{logical. Should the bias corrected version of misclassification
    error be computed?}
  \item{predictions}{logical, return a matrix of predictions. The ith column
contains predictions of the ith out-of-bootstrap sample and 'NA's
corresponding to the ith bootstrap sample.}
  \item{list.tindx}{list of numeric vectors, indicating which
    observations are included in each bootstrap sample.}
  \item{both.boot}{logical, return both (bootstrap and 632plus) estimations or
    only one of them.}
  \item{\dots}{additional arguments to \code{model}.}
}

\details{
 See \code{\link{errorest}}.
}
\keyword{misc}
