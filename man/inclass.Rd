\name{inclass}
\alias{inclass}
\alias{inclass.formula}
\alias{inclass.default}
\alias{inclass.flist}

\title{Indirect Classification}

\description{
A framework for the indirect classification approach.
}

\usage{
\method{inclass}{default}(object, ...)
\method{inclass}{formula}(formula, pFUN, data, subset, na.action, ...)
\method{inclass}{flist}(object, pFUN, data, subset, na.action, ...)
}

\arguments{
  \item{formula}{formula or object. A \code{formula} specified as \code{x+y+z~a+b+c}
is transformed into a list of formulas, where each vector \code{x, y, z} is modelled with formula
\code{a+b+c}. \code{object} is of class \code{\link{flist}}. }
  \item{pFUN}{function to be used to compute predictive models. Evaluated functions are: \code{rpart} - regression or classification trees (see
\code{\link[rpart]{rpart}}); \code{lm} - linear models (see \code{\link{lm}}); \code{bagging} - bagged regression or classification trees (see \code{\link{bagging}}).}
  \item{data}{data frame of explanatory and intermediate variables.}
  \item{subset}{optional vector specifying a subset of observations to be
used.}
  \item{na.action}{function which indicates what should happen when the data
contain \code{NA}s. Defaults to \code{getOption("na.action")}}
  \item{\dots}{additional arguments. For example \code{nbaggs} for \code{pFUN=bagging}, \ cr
          or \code{na.action=na.exclude} for
\code{pFUN = lm}.}
  \item{object}{ object of class \code{\link{flist}}.}
}

\details{
A given data set is subdivided into three classes of variables: those to be
used predicting the class (explanatory variables) those to be used defining
the class (intermediate variables) and the class membership variable itself
(response variable). Intermediate variables are modelled based on the explanatory variables.\cr 
Here, each specified intermediate variable is modelled separately 
following \code{pFUN} and a formula specified by \code{formula} 
or \code{flist}. A list of prediction models corresponding to each 
intermediate is returned. For a detailed description on indirect
classification see Hand et al. (2001).
}

\value{
 An object of class \code{inclass}, consisting of a list of prediction models for
 each specified intermediate variable. If \code{pFUN=bagging}, each element of the list (for each intermediate variable) consists of a list of \code{nbagg} models corresponding to each bootstrap sample.
}

\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics \& Data Analysis} \bold{36},
209--225.
}

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\seealso{\code{\link[rpart]{rpart}}, \code{\link{bagging}}, \code{\link{lm}}}

\examples{

data(smoking)
# Set three groups of variables:
# 1) explanatory variables are: TarY, NicY, COY, Sex, Age
# 2) intermediate variables are: TVPS, BPNL, COHB
# 3) response (resp) is defined by:

resp <- function(data){
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- resp(smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- cbind(smoking, response)

formula <- TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

inclass(formula, pFUN = lm, data = smoking)


}

\keyword{misc}
