\name{bagging}
\alias{bagging}
\alias{ipredbagg}
\alias{ipredbagg.factor}
\alias{ipredbagg.numeric}
\alias{ipredbagg.Surv}
\alias{ipredbagg.default}
\alias{bagging.data.frame}
\alias{bagging.default}
\title{Bagging Classification, Regression and Survival Trees }
\description{
  Bagging for classification, regression and survival trees.
}
\usage{
ipredbagg.factor(y, X=NULL, nbagg=25, control=
                 rpart.control(minsplit=2, cp=0, xval=0), 
                 comb=NULL, coob=FALSE, ns=length(y), keepX = TRUE, \dots)
ipredbagg.numeric(y, X=NULL, nbagg=25, control=rpart.control(xval=0), 
                  comb=NULL, coob=FALSE, ns=length(y), keepX = TRUE, \dots)
ipredbagg.Surv(y, X=NULL, nbagg=25, control=rpart.control(xval=0), 
               comb=NULL, coob=FALSE, ns=dim(y)[1], keepX = TRUE, \dots)
\method{bagging}{data.frame}(formula, data, subset, na.action=na.rpart, \dots)
}
\arguments{
  \item{y}{the response variable: either a factor vector of class labels
           (bagging classification trees), a vector of numerical values 
           (bagging regression trees) or an object of class 
           \code{\link[survival]{Surv}} (bagging survival trees).}
  \item{X}{a data frame of predictor variables.}
  \item{nbagg}{an integer giving the number of bootstrap replications. }
  \item{coob}{a logical indicating whether an out-of-bag estimate of the
              error rate (misclassification error, root mean squared error
              or Brier score) should be computed. 
              See \code{\link{predict.classbagg}} for
              details.}
  \item{control}{options that control details of the \code{rpart}
                algorithm, see \code{\link[rpart]{rpart.control}}. It is
                wise to set \code{xval = 0} in order to save computing 
                time. Note that the 
                default values depend on the class of \code{y}.}
  \item{comb}{a list of additional models for model combination, see below
for some examples. Note that argument \code{method} for double-bagging is no longer there, 
\code{comb} is much more flexible.}
  \item{ns}{number of sample to draw from the learning sample. By default,
            the usual bootstrap n out of n with replacement is performed. 
            If \code{ns} is smaller than \code{length(y)}, subagging
            (Buehlmann and Yu, 2002), i.e. sampling \code{ns} out of
            \code{length(y)} without replacement, is performed.}
  \item{keepX}{a logical indicating whether the data frame of predictors
    should be returned. Note that the computation of the 
    out-of-bag estimator requires  \code{keepX=TRUE}.}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} 
                 is the response variable and \code{rhs} a set of
                 predictors.}
  \item{data}{optional data frame containing the variables in the
              model formula.} 
  \item{subset}{optional vector specifying a subset of observations
                to be used.}
  \item{na.action}{function which indicates what should happen when
                   the data contain \code{NA}s.  Defaults to
                   \code{\link[rpart]{na.rpart}}.}
  \item{...}{additional parameters passed to \code{ipredbagg} or 
\code{\link[rpart]{rpart}}, respectively.}
}

\details{
  Bagging for classification and regression trees were suggested by
Breiman (1996a, 1998) in order to stabilise trees. 

The trees in this function are computed using the implementation in the 
\code{\link[rpart]{rpart}} package. The generic function \code{ipredbagg}
implements methods for different responses. If \code{y} is a factor,
classification trees are constructed. For numerical vectors
\code{y}, regression trees are aggregated and if \code{y} is a survival 
object, bagging survival trees (Hothorn et al, 2002) is performed. 
The function \code{bagging} offers a formula based interface to
\code{ipredbagg}.

\code{nbagg} bootstrap samples are drawn and a tree is constructed 
for each of them. There is no general rule when to stop the tree 
growing. The size of the
trees can be controlled by \code{control} argument 
or \code{\link{prune.classbagg}}. By
default, classification trees are as large as possible whereas regression
trees and survival trees are build with the standard options of
\code{\link[rpart]{rpart.control}}. If \code{nbagg=1}, one single tree is
computed for the whole learning sample without bootstrapping.

If \code{coob} is TRUE, the out-of-bag sample (Breiman,
1996b) is used to estimate the prediction error 
corresponding to \code{class(y)}. Alternatively, the out-of-bag sample can
be used for model combination, an ouf-of-bag error rate estimator is not 
available in this case. Double-bagging (Hothorn and Lausen,
2002a) computes a LDA on the out-of-bag sample and uses the discriminant
variables as additional predictors for the classification trees. \code{comb}
is an optional list of lists with two elements \code{model} and \code{predict}. 
\code{model} is a function with arguments \code{formula} and \code{data}. 
\code{predict} is a function with arguments \code{object, newdata} only. If
the estimation of the covariance matrix in \code{\link{lda}} fails due to a
limited out-of-bag sample size, one can use \code{\link{slda}} instead.
See the example section for an example of double-bagging. The methodology is
not limited to a combination with LDA: bundling (Hothorn and Lausen, 2002b) 
can be used with arbritrary classifiers.

}

\value{
  The class of the object returned depends on \code{class(y)}:
\code{classbagg, regbagg} and \code{survbagg}. Each is a list with elements
\item{y}{the vector of responses.}
\item{X}{the data frame of predictors.}
\item{mtrees}{multiple trees: a list of length \code{nbagg} containing the
trees (and possibly additional objects) for each bootstrap sample.}
\item{OOB}{logical whether the out-of-bag estimate should be computed.}
\item{err}{if \code{OOB=TRUE}, the out-of-bag estimate of
misclassification or root mean squared error or the Brier score for censored
data.}
\item{comb}{logical whether a combination of method was requested.}

 For each class methods for the generics \code{\link[rpart]{prune}}, 
\code{\link{print}}, \code{\link{summary}} and \code{\link{predict}} are
available for inspection of the results and prediction, for example:
\code{\link{print.classbagg}}, \code{\link{summary.classbagg}}, 
\code{\link{predict.classbagg}}  and \code{\link{prune.classbagg}} for
classification problems.

}

\references{ 

Leo Breiman (1996a), Bagging Predictors. \emph{Machine Learning}
\bold{24}(2), 123--140.

Leo Breiman (1996b), Out-Of-Bag Estimation. \emph{Technical Report}
\url{ftp://ftp.stat.berkeley.edu/pub/users/breiman/OOBestimation.ps.Z}.

Leo Breiman (1998), Arcing Classifiers. \emph{The Annals of Statistics}
\bold{26}(3), 801--824.

Peter Buehlmann and Bin Yu (2002), Analyzing Bagging. \emph{The Annals of
Statistics} \bold{30}(4), 927--961.

Torsten Hothorn and Berthold Lausen (2002a), Double-Bagging: Combining
classifiers by bootstrap aggregation. \emph{Pattern Recognition}, to appear.
Preprint available from 
\url{http://www.mathpreprints.com/math/Preprint/hothorn/20020227.2/1}.

Torsten Hothorn and Berthold Lausen (2002b), Bundling Classifiers by Bagging
Trees. \emph{submitted}.
Preprint available from 
\url{http://www.mathpreprints.com/math/Preprint/blausen/20021016/1}.

Torsten Hothorn, Berthold Lausen, Axel Benner and Martin
Radespiel-Troeger (2002), Bagging Survival Trees. \emph{submitted}.
Preprint available from
\url{http://www.mathpreprints.com/math/Preprint/blausen/20020518/2}.

}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

# Classification: Breast Cancer data

data(BreastCancer)

# Test set error bagging (nbagg = 50): 3.7\% (Breiman, 1998, Table 5)

mod <- bagging(Class ~ Cl.thickness + Cell.size
                + Cell.shape + Marg.adhesion   
                + Epith.c.size + Bare.nuclei   
                + Bl.cromatin + Normal.nucleoli
                + Mitoses, data=BreastCancer, coob=TRUE)
print(mod)

# Test set error bagging (nbagg=50): 7.9\% (Breiman, 1996a, Table 2)

data(Ionosphere)
Ionosphere$V2 <- NULL # constant within groups

bagging(Class ~ ., data=Ionosphere, coob=TRUE)

# Double-Bagging: combine LDA and classification trees

# predict returns the linear discriminant values, i.e. linear combinations
# of the original predictors

comb.lda <- list(list(model=lda, predict=function(obj, newdata)
                                 predict.lda(obj, newdata)$x))

# Note: out-of-bag estimator is not available in this situation, use
# errorest

mod <- bagging(Class ~ ., data=Ionosphere, comb=comb.lda) 

predict(mod, Ionosphere[1:10,])

# Regression:


data(BostonHousing)

# Test set error (nbagg=25, trees pruned): 3.41 (Breiman, 1996a, Table 8)

mod <- bagging(medv ~ ., data=BostonHousing, coob=TRUE)
print(mod)

learn <- as.data.frame(mlbench.friedman1(200))

# Test set error (nbagg=25, trees pruned): 2.47 (Breiman, 1996a, Table 8)

mod <- bagging(y ~ ., data=learn, coob=TRUE)
print(mod)

# Survival data

# Brier score for censored data estimated by 
# 10 times 10-fold cross-validation: 0.2 (Hothorn et al,
# 2002)

data(DLBCL)
mod <- bagging(Surv(time,cens) ~ MGEc.1 + MGEc.2 + MGEc.3 + MGEc.4 + MGEc.5 +
                                 MGEc.6 + MGEc.7 + MGEc.8 + MGEc.9 +
                                 MGEc.10 + IPI, data=DLBCL, coob=TRUE)

print(mod)


}
\keyword{tree}
