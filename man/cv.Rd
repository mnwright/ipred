\name{cv}
\alias{cv}
\alias{cv.default}
\alias{cv.factor}
\alias{cv.numeric}
\alias{cv.Surv}
\title{Cross-validated Error Rate Estimators.}
\description{
  Those functions are low-level functions used by \code{\link{errorest}} and
are normally not called by users.
}
\usage{
\method{cv}{factor}(y, X, model, predict, k=10, random=TRUE, strat=TRUE,
            predictions=FALSE, iformula=NULL, \dots) 
}
\arguments{
  \item{y}{response variable, either of class \code{factor}
(classification), \code{numeric} (regression) or \code{Surv} (survival).}
  \item{X}{data frame of predictors.}
  \item{model}{a function implementing the predictive model to be
evaluated.}
  \item{predict}{a function with arguments \code{object} and \code{newdata}
only which predicts the status of the observations in \code{newdata} based
on the fitted model in \code{object}.}
  \item{k}{k-fold cross-validation.}
  \item{random}{logical, indicates whether a random order or the given
order of the data should be used for sample splitting or not.}
  \item{strat}{logical, stratified sampling or not.}
  \item{predictions}{logical, return the prediction of each observation.}
  \item{iformula}{formula for indirect classification.}
  \item{\dots}{additional arguments to \code{model}.}
}
\details{
 See \code{\link{errorest}}.
}
\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{misc}
