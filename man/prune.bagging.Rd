\name{prune.bagging}
\alias{prune.bagging}
\title{ Prune Bagged Trees  }
\description{
  Prune each of the trees in a \code{bagging} object. 
}
\usage{
\method{prune}{bagging}(tree, cp=0.01,...)
}
\arguments{
  \item{tree}{ object of class \code{bagging} (calling this \code{tree}
is needed by the generic method \code{prune} in package \code{rpart}).}
  \item{cp}{complexity parameter, see \code{\link[rpart]{prune.rpart}}.}
  \item{...}{additional arguments to \code{\link[rpart]{prune.rpart}}.}
}
\details{
  By default, \code{\link{bagging}} grows trees of maximal size. One may
want to prune each tree, however, it is not clear whether or not 
this may decrease prediction error. 
}
\value{
  An object of class bagging with the pruned trees is returned. 
}
\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))

mt <- bagging(y, X)
pt <- prune(mt, cp=0.01)

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))

cls <- predict(mt, X)

cat("error: ", sum(y != cls)/length(cls), "\n")

cls <- predict(pt, X)

cat("error: ", sum(y != cls)/length(cls), "\n")


}
\keyword{tree}
