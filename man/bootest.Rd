\name{bootest}
\alias{bootest}
\alias{bootest.default}
\alias{bootest.factor}
\alias{bootest.numeric}
\alias{bootest.Surv}
\title{Bootstrap Error Rate Estimators}
\description{
  Those functions are low-level functions used by \code{\link{errorest}} and
are normally not called by users.
}
\usage{
\method{bootest}{factor}(y, X, model, predict, nboot=25, bc632plus=FALSE,
iformula=NULL, \dots)
}
\arguments{
  \item{y}{the response variable, either of class \code{factor}
(classification), \code{numeric} (regression) or \code{Surv} (survival).}
  \item{X}{the data frame of predictors.}
  \item{model}{a function implementing the predictive model to be
evaluated.}
  \item{predict}{a function with arguments \code{object} and \code{newdata}
only which predicts the status of the observations in \code{newdata} based
on the fitted model in \code{object}.}
  \item{nboot}{number of bootstrap replications to be used.}
  \item{bc632plus}{logical. Should the bias corrected version of misclassification
error be computed?}
 \item{iformula}{formula for indirect classification.}

  \item{\dots}{additional arguments to \code{model}.}
}

\details{
 See \code{\link{errorest}}.
}
\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{misc}
