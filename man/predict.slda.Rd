\name{predict.slda}
\alias{predict.slda}
\title{ Predictions from Stabilised Linear Discriminant Analysis }
\description{
  Predict the class of a new observation based on stabilised LDA.  
}
\usage{
\method{predict}{slda}(object, newdata, ...)
}
\arguments{
  \item{object}{object of class \code{slda}.}
  \item{newdata}{a data frame of new observations. }
  \item{...}{additional arguments passed to
\code{\link[MASS]{predict.lda}}.}
}
\details{
  This function is a method for the generic function \code{\link{predict}}
for class \code{slda}. For the details see \code{\link[MASS]{predict.lda}}.

}
\value{
A list with components
\item{class}{the predicted class (a factor).}
\item{posterior}{posterior probabilities for the classes.}
\item{x}{the scores of test cases.}
}


\keyword{multivariate}
