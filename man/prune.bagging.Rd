\name{prune.classbagg}
\alias{prune.classbagg}
\alias{prune.regbagg}
\alias{prune.survbagg}
\title{ Pruning for Bagging  }
\description{
  Prune each of the trees returned by \code{\link{bagging}}. 
}
\usage{
\method{prune}{classbagg}(tree, cp=0.01,...)
}
\arguments{
  \item{tree}{ an object returned by \code{\link{bagging}} 
              (calling this \code{tree} is needed by the generic function 
              \code{prune} in package \code{rpart}).}
  \item{cp}{complexity parameter, see \code{\link[rpart]{prune.rpart}}.}
  \item{...}{additional arguments to \code{\link[rpart]{prune.rpart}}.}
}
\details{
  By default, \code{\link{bagging}} grows classification 
  trees of maximal size. One may want to prune each tree, however, 
  it is not clear whether or not this may decrease prediction error. 
}
\value{
  An object of the same class as \code{tree} with the trees pruned. 
}
\examples{

data("Glass", package = "mlbench")

mod <- bagging(Type ~ ., data=Glass, nbagg=10, coob=TRUE)
pmod <- prune(mod)
print(pmod)


}
\keyword{tree}
