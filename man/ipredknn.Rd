\name{ipredknn}
\alias{ipredknn}
\title{ k-Nearest Neighbour Classification }
\description{
  $k$-nearest neighbour classification with an interface compatible to 
\code{\link{bagging}} and \code{\link{errorest}}.
}
\usage{
ipredknn(formula, data, subset, na.action, k=5, \dots)
}
\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} 
                 is the response variable and \code{rhs} a set of
                 predictors.}
  \item{data}{optional data frame containing the variables in the
              model formula.} 
  \item{subset}{optional vector specifying a subset of observations
                to be used.}
  \item{na.action}{function which indicates what should happen when
                   the data contain \code{NA}s.}
  \item{k}{number of neighbours considered, defaults to 5.}
  \item{...}{additional parameters.}
}

\details{
  This is a wrapper to \code{\link[class]{knn}} in order to be able to 
  use k-NN in \code{\link{bagging}} and \code{\link{errorest}}. 
}

\value{
  An object of class \code{ipredknn}. See \code{\link{predict.ipredknn}}.
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

learn <- as.data.frame(mlbench.twonorm(300))

mypredict.knn <- function(object, newdata) 
                   predict.ipredknn(object, newdata, type="class")

errorest(classes ~., data=learn, model=ipredknn, 
         predict=mypredict.knn)


}
\keyword{multivariate}
