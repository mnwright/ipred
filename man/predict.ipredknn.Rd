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
  \item{type}{return either the predicted class or the 
              the proportion of the votes for the winning
              class.}
  \item{...}{additional arguments.}
}
\details{
  This function is a method for the generic function \code{\link{predict}}
for class \code{ipredknn}. For the details see \code{\link[class]{knn}}.

}
\value{
Either the predicted class or the
the proportion of the votes for the winning class.
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\keyword{multivariate}
