\name{cv}
\alias{cv}
\alias{cv.default}
\alias{cv.factor}
\alias{cv.numeric}
\alias{cv.integer}
\alias{cv.Surv}
\title{Cross-validated Error Rate Estimators.}
\description{
  Those functions are low-level functions used by \code{\link{errorest}} and
are normally not called by users.
}
\usage{
\method{cv}{factor}(y, formula, data, model, predict, k=10, random=TRUE, 
            strat=FALSE,
            predictions=NULL, getmodels=NULL, list.tindx = NULL, \dots) 
}
\arguments{
  \item{y}{response variable, either of class \code{factor}
(classification), \code{numeric} (regression) or \code{Surv} (survival).}
  \item{formula}{a formula object.}
  \item{data}{data frame of predictors and response described in \code{formula}.}
  \item{model}{a function implementing the predictive model to be
evaluated. The function \code{model} can either return an
               object representing a fitted model or a function with
               argument \code{newdata} which returns predicted values. In
               this case, the \code{predict} argument to \code{errorest} is
               ignored.}
  \item{predict}{a function with arguments \code{object} and \code{newdata}
only which predicts the status of the observations in \code{newdata} based
on the fitted model in \code{object}.}
  \item{k}{k-fold cross-validation.}
  \item{random}{logical, indicates whether a random order or the given
order of the data should be used for sample splitting or not, defaults to
\code{TRUE}.}
  \item{strat}{logical, stratified sampling or not, defaults to \code{FALSE}.}
  \item{predictions}{logical, return the prediction of each observation.}
  \item{getmodels}{logical, return a list of models for each fold.}
  \item{list.tindx}{list of numeric vectors, indicating which
    observations are included in each cross-validation sample.}
  \item{\dots}{additional arguments to \code{model}.}
}
\details{
 See \code{\link{errorest}}.
}
\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{misc}
