\name{inclass}
\alias{inclass}
\alias{inclass.default}
\alias{inclass.data.frame}

\title{Indirect Classification}

\description{
A framework for the indirect classification approach.
}

\usage{
\method{inclass}{data.frame}(formula, data, pFUN = NULL, cFUN = NULL, ...)
}
\arguments{
  \item{formula}{formula. A \code{formula} specified as
    \code{y~w1+w2+w3~x1+x2+x3} models each intermediate variable
    \code{w1, w2, w3} by \code{wi~x1+x2+x3} and the response by
    \code{y~w1+w2+w3} if no other formulas are given in \code{pFUN} or \code{cFUN}.}
  \item{data}{data frame of explanatory, intermediate and response variables.}
  \item{pFUN}{list of lists, which describe models for the intermediate variables, see below for details.}
  \item{cFUN}{either a function or a list which describes the model for the
response variable. The function has the argument \code{newdata} only.}
  \item{\dots}{additional arguments, passed to model fitting of the
    response variable.} 
}

\details{
A given data set is subdivided into three types of variables: those to be
used predicting the class (explanatory variables) those to be used defining
the class (intermediate variables) and the class membership variable itself
(response variable). Intermediate variables are modelled based on the
explanatory variables, the class membership variable is defined on the
intermediate variables.\cr

Each specified intermediate variable is modelled separately 
following \code{pFUN} and a formula specified by \code{formula}.
\code{pFUN} is a list of lists, the maximum length of
\code{pFUN} is the number of intermediate variables. Each element of
\code{pFUN} is a list with elements:\cr
  \code{model} -  a function with arguments \code{formula} and
  \code{data}; \cr
  \code{predict} - an optional function with arguments \code{object, newdata} only, 
  if \code{predict} is not specified, the predict method of \code{model}
  is used; \cr
  \code{formula} - specifies the formula for the corresponding
  \code{model} (optional),
    the formula described in \code{y~w1+w2+w3~x1+x2+x3} is used if no other is
  specified.
\cr
  
The response is classified following \code{cFUN}, which is either a fixed
function or a list as described below. The determined function \code{cFUN} assigns the intermediate (and
explanatory) variables to a certain class membership, the list
\code{cFUN} has the elements \code{formula, model, predict} and
\code{training.set}. The elements \code{formula, model, predict} are
structured as described by \code{pFUN}, the described model is
trained on the original (intermediate variables) if \code{training.set="original"}
or if \code{training.set = NULL}, on the fitted values if
\code{training.set = "fitted"} or on observations not included in a
specified subset if \code{training.set = "subset"}.
\cr 

A list of prediction models corresponding to each 
intermediate variable, a predictive function for the response, a list of
specifications for the intermediate and for the response are returned. \cr
For a detailed description on indirect
classification see Hand et al. (2001).
}
\value{
 An object of class \code{inclass}, consisting of a list of 
  \item{model.intermediate}{list of fitted models for each intermediate
variable.}
  \item{model.response}{predictive model for the response variable.}
  \item{para.intermediate}{list, where each element is again a list and specifies
the model for each intermediate variable.}
  \item{para.response}{a list which specifies the model for response variable.}
}
\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics & Data Analysis} \bold{36},
209--225.

Andrea Peters, Berthold Lausen, Georg Michelson and Olaf Gefeller (2003),
Diagnosis of glaucoma by indirect classifiers.
\emph{Methods of Information in Medicine} \bold{1}, 99-103.
}

\seealso{\code{\link{bagging}}, \code{\link{inclass}}}

\examples{
data("Smoking", package = "ipred")
# Set three groups of variables:
# 1) explanatory variables are: TarY, NicY, COY, Sex, Age
# 2) intermediate variables are: TVPS, BPNL, COHB
# 3) response (resp) is defined by:

classify <- function(data){
  data <- data[,c("TVPS", "BPNL", "COHB")]
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- classify(Smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- data.frame(Smoking, response)

formula <- response~TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

inclass(formula, data = smoking, pFUN = list(list(model = lm, predict =
mypredict.lm)), cFUN = classify)

}

\keyword{misc}
