\name{flist}
\alias{flist}

\title{Create an Object of Class flist }
\description{
 An object of class \code{flist} is created.
}
\usage{
flist(formula, ...)
}

\arguments{
  \item{formula}{formula or list of formulas.}
  \item{\dots}{additional formulas. }
}
\value{
 \code{flist} is a list of formulas, where each element has one response variable. 
  An object of class \code{flist} is returned, see \code{\link{inclass}}.
}

\note{
 If \code{formula} is a list, only one response variable per element can be
specified. The additional formulas in \code{\dots} must also have only one
response variable. A \code{formula} specified as \code{x+y+z~a+b+c} is
transformed into a list with one response variable per element and added to
the list of the formulas specified in \code{dots}.
}

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\seealso{\code{\link{inclass.flist}}}

\keyword{misc}
