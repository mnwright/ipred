\name{slda}
\alias{slda}
\alias{slda.default}
\alias{slda.formula}
\alias{slda.factor}
\title{ Stabilised Linear Discriminant Analysis }
\description{
  Linear discriminant analysis based on left-spherically 
  distributed linear scores. }
}
\usage{
\method{slda}{formula}(formula, data, subset, na.action=na.rpart, \dots)
\method{slda}{factor}(y, X, q=NULL, \dots)
}
\arguments{
  \item{y}{the response variable: a factor vector of class labels.}
  \item{X}{a data frame of predictor variables.}
  \item{q}{the number of positive eigenvalues the scores are derived from,
           see below.}
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
  \item{...}{additional parameters passed to \code{\link[MASS]{lda}}.}
}

\details{
  This function implements the LDA for \eqn{q}-dimensional linear scores of
the original \eqn{p} predictors derived from the \eqn{PC_q} rule by Laeuter
et al. (1998). Based on the product sum matrix 
\deqn{W = (X - \bar{X})^\top(X - \bar{X})}
the eigenvalue problem \eqn{WD = diag(W)DL} is solved. The first \eqn{q}
columns \eqn{D_q} of \eqn{D} are used as a weight matrix for the 
original \eqn{p} predictors: \eqn{XD_q}. By default, \eqn{q} is the number
of eigenvalues greater one. The \eqn{q}-dimensional linear scores are
left-spherically distributed and are used as predictors for a classical 
LDA. 

This form of reduction of the dimensionality was 
developed for discriminant analysis problems by Laeuter (1992) and was used
for multivariate tests by Laeuter et al. (1998), Kropf (2000) gives an
overview. For details on left-spherically distributions see Fang and 
Zhang (1990).  

}

\value{
  An object of class \code{slda}, a list with components
  \item{scores}{the weight matrix.}
  \item{mylda}{an object of class \code{lda}.}
}

\seealso{
  \code{\link{predict.slda}}
}

\references{ 

Fang Kai-Tai and Zhang Yao-Ting (1990), \emph{Generalized Multivariate
Analysis}, Springer, Berlin.

Siegfried Kropf (2000), \emph{Hochdimensionale multivariate Verfahren in der
medizinischen Statistik}, Shaker Verlag, Aachen (in german).

Juergen Laeuter (1992), \emph{Stabile multivariate Verfahren},
Akademie Verlag, Berlin (in german).

Juergen Laeuter, Ekkehard Glimm and Siegfried Kropf (1998), Multivariate
Tests Based on Left-Spherically Distributed Linear Scores. \emph{The Annals
of Statistics}, \bold{26}(5) 1972--1988. 



}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

learn <- as.data.frame(mlbench.twonorm(100))
test <- as.data.frame(mlbench.twonorm(1000))

mlda <- lda(classes ~ ., data=learn)
mslda <- slda(classes ~ ., data=learn)

print(mean(predict(mlda, newdata=test)$class != test$classes))
print(mean(predict(mslda, newdata=test)$class != test$classes))

}
\keyword{multivariate}
