\name{summary.classbagg}
\alias{summary.classbagg}
\alias{summary.regbagg}
\alias{summary.survbagg}
\alias{print.summary.bagging}
\title{Summarising Bagging}
\description{
  \code{summary} method for objects returned by \code{\link{bagging}}.
}
\usage{
\method{summary}{classbagg}(object, \dots)
}
\arguments{
  \item{object}{object returned by \code{\link{bagging}}.}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
\details{
  A representation of all trees in the object is printed.
}
\value{
  none
}
\keyword{tree}
