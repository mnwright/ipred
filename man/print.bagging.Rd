\name{print.bagging}
\alias{print.bagging}
\title{Print Method for Bagging Object}
\description{
  Print object of class \code{"bagging"} in nice layout.
}
\usage{
\method{print}{bagging}(x, digits, \dots)
}
\arguments{
  \item{x}{object of class \code{bagging}.}
  \item{digits}{how many digits should be printed.}
  \item{\dots}{further arguments to be passed to or from methods.}  
}
\details{
  A \code{bagging} object is a named list of numbers and
  character strings, supplemented with \code{method} and \code{note}
  elements.  The \code{method} is displayed as a title, the \code{note}
  as a footnote and the remaining elements are given in an aligned
  `name = value' format.
}
\value{
  none
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\keyword{tree}
