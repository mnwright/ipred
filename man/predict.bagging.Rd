\name{predict.bagging}
\alias{predict.bagging}
\title{ Predictions from a Bagging Object }
\description{
  Predict the outcome of a new observation using bagged trees.  
}
\usage{
predict.bagging(object, newdata=NULL,...)
}
\arguments{
  \item{object}{ object of class \code{\link{bagging}}.}
  \item{newdata}{ new observations. }
  \item{...}{additional arguments.}
}
\details{
  The predicted class is derived from majority voting over the classes
predicted by all classification trees. In the regression framework, the mean
of the predictions is returned. 

}
\value{
  The predicted class or endpoint in regression of all observations are returned.
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\references{ 

Leo Breiman (1996), Bagging Predictors. \emph{Machine Learning}
\bold{24}(2), 123--140.

}

\examples{

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))

mt <- bagging(y, X)

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))

cls <- predict(mt, X)

cat("error: ", sum(y != cls)/length(cls), "\n")


}
\keyword{tree}
