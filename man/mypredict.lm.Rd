\name{mypredict.lm}
\alias{mypredict.lm}
\title{Predictions based on Linear Models}
\description{
Function to predict a vector of full length (number of observations), where predictions according to missing
explanatory values are replaced by \code{NA}.
}

\usage{
mypredict.lm(object, newdata)
}

\arguments{
  \item{object}{an object of class \code{lm}.}
  \item{newdata}{matrix or data frame to be predicted according to \code{object}.}
}

\value{
Vector of predicted values.
}

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\note{\code{predict.lm} delivers a vector of reduced length, i.e. rows where
explanatory variables are missing are omitted. The full length of the
predicted observation vector is necessary in the indirect classification
approach (\code{\link{predict.inclass}}).}

\keyword{misc}
