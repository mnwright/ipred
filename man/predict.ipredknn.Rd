\name{predict.ipredknn}
\alias{predict.ipredknn}
\title{ Predictions from k-Nearest Neighbors }
\description{
  Predict the class of a new observation based on k-NN. 
}
\usage{
\method{predict}{ipredknn}(object, newdata, type=c("prob", "class"), ...)
}
\arguments{
  \item{object}{object of class \code{ipredknn}.}
  \item{newdata}{a data frame of new observations. }
  \item{type}{return either class probability estimates or the predicted
class itself.}
  \item{...}{additional arguments.}
}
\details{
  This function is a method for the generic function \code{\link{predict}}
for class \code{ipredknn}. For the details see \code{\link[class]{knn}}.

}
\value{
  Either a matrix of predicted class probabilities or the predicted classes
itself.
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\keyword{multivariate}
