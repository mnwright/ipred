\name{bagging}
\alias{bagging}
\alias{bagging.formula}
\alias{bagging.default}
\title{Bagging Classification and Regression Trees }
\description{
  Bootstrap aggregated classification and regression trees.
}
\usage{
\method{bagging}{default}(y, X=NULL, nbagg=25, method=c("standard","double"),
        coob=TRUE, control= rpart.control(minsize=2, cp=0), \dots)
\method{bagging}{formula}(formula, data, subset, na.action=na.rpart, \dots)
}
\arguments{
  \item{y}{ vector of responses: either numerical (regression) or factors (classification). }
  \item{X}{ data frame of predictors. }
  \item{nbagg}{ number of bootstrap replications. }
  \item{method}{ \code{standard} for Bagging and \code{double} for
Double-Bagging. }
  \item{coob}{logical. Compute an out-of-bag estimate of the
misclassification or mean-squared error.}
  \item{control}{options that control details of the \code{rpart}
algorithm, see \code{\link[rpart]{rpart.control}}.}
  \item{formula}{ formula describing the model: \code{y ~ x + w + z}, where \code{y} is the response 
  and \code{x,w,z} are predictors, see \code{\link{lm}} for details.}
  \item{data}{ optional data frame containing the variables in the
    model formula.}
  \item{subset}{ optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{ function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{\link[rpart]{na.rpart}}.}
  \item{...}{additional parameters to methods (e.g. \code{\link[rpart]{rpart}}).}
}

\details{
  Bootstrap aggregated classification and regression trees were suggested by
Breiman (1996, 1998) in order to stabilise trees. This function
is based on trees computed by \code{\link[rpart]{rpart}}. If \code{y} is a
factor, classification trees are constructed, regression trees otherwise. 
\code{nbagg} bootstrap samples are drawn and a tree is constructed 
for each of them. If \code{coob} is TRUE, the out-of-bag sample is 
used to estimate the prediction error. Double-Bagging (Hothorn and Lausen,
2002) computes a LDA on the out-of-bag sample and uses the discriminant
variables as additional predictors for the classification trees. Therefore,
an out-of-bag estimate of misclassification error is not available for
\code{method="double"}.

\code{\link{print.bagging}} and
\code{\link{summary.bagging}} are available for the inspection of the
results as well as \code{\link{predict.bagging}} for prediction.
Additionally, the function \code{\link{prune.bagging}} can be used to prune
each of the \code{nbagg} trees. By default, the trees are not pruned and the
tree growing is not stopped until the nodes are pure.
}

\value{
  An object of class \code{bagging}: a list containing the following objects
  \item{mt}{list of length \code{nbagg} containing \code{\link[rpart]{rpart}}
trees.}
  \item{oob}{out-of-bag predictions for each observation.}
  \item{err}{out-of-bag error estimate.}
  \item{nbagg}{number of bootstrap samples and trees used.}
  \item{method}{method used.}
  \item{ldasc}{discriminant functions of LDA (for Double-Bagging only).}
}

\references{ 

Leo Breiman (1996), Bagging Predictors. \emph{Machine Learning}
\bold{24}(2), 123--140.

Leo Breiman (1998), Arcing Classifiers. \emph{The Annals of Statistics}
\bold{26}(3), 801--824.

Torsten Hothorn and Berthold Lausen (2002), Double-Bagging: Combining
classifiers by bootstrap aggregation. \emph{submitted}, 
preprint available under 
\url{http://www.mathpreprints.com/math/Preprint/hothorn/20020227.2/1}.

}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))
learn <- cbind(y, X)

mt <- bagging(y ~., data=learn, coob=TRUE)
mt

dmt <- bagging(y ~., data=learn, coob=FALSE, method="double")
dmt

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))

cls <- predict(mt, newdata=X)

cat("Misclass error est: ", mean(y != cls), "\n")
cat("Misclass error oob: ", mt$err, "\n")

cls <- predict(dmt, newdata=X)

cat("Misclass error est double-bagging: ", mean(y != cls), "\n")

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- apply(X, 1, mean) + rnorm(nrow(X))

learn <- cbind(y, X)

mt <- bagging(y ~., data=learn, coob=TRUE)
mt
 
X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- apply(X, 1, mean) + rnorm(nrow(X))

haty <- predict(mt, newdata=X)

cat("MSE error: ", mean((haty - y)^2) , "\n")

data(BreastCancer)
BreastCancer$Id <- NULL

# Test set error bagging (nbagg = 50): 3.7\% (Breiman, 1998, Table 5)

mod <- bagging(Class ~ Cl.thickness + Cell.size 
                + Cell.shape + Marg.adhesion
                + Epith.c.size + Bare.nuclei
                + Bl.cromatin + Normal.nucleoli
                + Mitoses, data=BreastCancer, coob=TRUE)
mod
# Do not cross-validate for performance reasons ;-)
print(mean(BreastCancer$Class != predict(mod, newdata=BreastCancer)))

mod <- bagging(Class ~ Cl.thickness + Cell.size 
                + Cell.shape + Marg.adhesion
                + Epith.c.size + Bare.nuclei
                + Bl.cromatin + Normal.nucleoli
                + Mitoses, data=BreastCancer, method="double")
mod
# Do not cross-validate for performance reasons ;-)
print(mean(BreastCancer$Class != predict(mod, newdata=BreastCancer)))

}
\keyword{tree}
