\name{errorest}
\alias{errorest}

\title{Estimators for the Prediction Error}
\description{
Resampling based estimates of prediction error (misclassification or 
mean squared error). 
}
\usage{
errorest(formula, data, subset, na.action, model=NULL, predict=NULL,
         iclass=NULL, estimator=c("cv", "boot", "632plus"), 
         est.para=list(k = 10, nboot = 25), ...)
}

\arguments{
  \item{formula}{formula. Either describing the model of explanatory and 
                 response variables in the usual way (see \code{\link{lm}}) 
                 or the model between explanatory and intermediate variables
                 in the framework of indirect classification, 
                 see \code{\link{inclass}}.}
  \item{data}{data frame containing the variables in the model formula 
              and additionally the class membership variable 
              if \code{model = inclass}. \code{data} is required for
              indirect classification, otherwise \code{formula} is evaluated
              in the calling environment.}
  \item{subset}{optional vector, specifying a subset of observations to 
                be used.}
  \item{na.action}{function. Indicates what should happen when the data 
                   contain \code{NA}s.}
  \item{model}{function. Modelling technique whose error rate is to be 
               estimated. The parameter \code{na.action} is evaluated in 
               the modelling process.}
  \item{predict}{function. Prediction method to be used. The vector of 
                 predicted values must have the same length as the the 
                 number of to-be-predicted observations. Predictions 
                 corresponding to missing data must be replaced by \code{NA}.
                 Additionally, \code{predict} has to return predicted values 
                 comparable to the responses (that is: factors for 
                 classification problems). See the example on how to make
                 this sure for any predictor.}
  \item{iclass}{character. Specifying the class membership variable 
                (response) in \code{data} in the framework of indirect
                classification.}
  \item{estimator}{estimator of the misclassification error: 
                   \code{cv} cross-validation, \code{boot} bootstrap or 
                   \code{632plus} bias corrected bootstrap (classification
                   only). }
  \item{est.para}{a list of additional parameters for the estimator. \code{k} 
                  for k-fold cross-validations or \code{nboot} for the number 
                  of bootstrap replications.}
  \item{\dots}{additional parameters to \code{model}.}
}
\details{
  The prediction error for classification and regression models using 
cross-validation or the bootstrap can be computed by \code{errorest}. Any
model can be specified as long as this is a function with arguments
\code{model(formula, data, subset, na.action, ...)}. If \code{model} is
generic and a \code{model.predict(object, newdata, ...)} is available,
\code{predict} does not need to be specified. However, \code{predict} has to
return predicted values directly comparable to the responses. See the
examples below. 

k-fold cross-validation and the usual bootstrap estimator with
\code{est.para$nboot} bootstrap replications can be computed for
classification and regression problems. The bias corrected .632+ bootstrap
by Efron and Tibshirani (1997) is available for classification problems
only. 

\code{\link{print.errorest}} is available for inspection of the results.

}

\value{
 An object of class \code{errorest}, i.e. a list with arguments:
  \item{err}{estimated misclassification error for a nominal or the square
root of the estimated mean squared error for a continuous response.}
  \item{estimator}{kind of estimator used.}
  \item{para}{additional parameters for the estimator.}
  \item{data.name}{names of the variables used.}
  \item{class}{logical. \code{TRUE} for classification problems.}
  \item{sd}{jackknife estimate of internal standard deviation of \code{err}
(if \code{estimator = "boot"}).}
}

\references{
Bradley Efron and Robert Tibshirani (1997),
Improvements on Cross-Validation: The .632+ Bootstrap Estimator.
\emph{Journal of the American Statistical Association} \bold{92}(438),
548--560.

Brian D. Ripley (1996), \emph{Pattern Recognition and Neural Networks},
Cambridge: Cambridge University Press.

David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics \& Data Analysis} \bold{36},
209--225.

}
\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de> and \cr
  Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>}
\examples{

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))
learn <- cbind(y, X)

mypredict.lda <- function(object, newdata)
  predict(object, newdata = newdata)$class

errorest(y ~ ., data= learn, model=lda, 
         estimator = "cv", predict= mypredict.lda)

# n-fold cv = leave-one-out.

errorest(y ~ ., data= learn, model=lda, 
         estimator = "cv", est.para=list(k = nrow(learn)), 
         predict= mypredict.lda)

errorest(y ~ ., data= learn, model=lda, 
         estimator = "boot", predict= mypredict.lda)

errorest(y ~ ., data= learn, model=lda, 
         estimator = "632plus", predict= mypredict.lda)

attach(learn)
errorest(y ~ V1 + V2 + V3, model=lda, estimator = "cv",
         predict= mypredict.lda)
detach(learn)


mypredict.rpart <- function(object, newdata)
  predict(object, newdata = newdata,type="class")

errorest(y ~ ., data= learn, model=rpart, estimator = "cv",
         predict=mypredict.rpart)

errorest(y ~ ., data= learn, model=rpart, estimator = "boot",
predict=mypredict.rpart)

errorest(y ~ ., data= learn, model=rpart, estimator = "632plus",
predict=mypredict.rpart)

errorest(y ~ ., data= learn, model=bagging, estimator = "cv",
nbagg=10)

data(Glass)

# LDA has cross-validated misclassification error of 
# 38\% (Ripley, 1996, page 98)


# Pruned trees about 32\% (Ripley, 1996, page 230)

pruneit <- function(formula, ...)
  prune(rpart(formula, ...), cp =0.01)

errorest(Type ~ ., data=Glass, model=pruneit, estimator= "cv",
predict=mypredict.rpart)

data(smoking)
# Set three groups of variables:
# 1) explanatory variables are: TarY, NicY, COY, Sex, Age
# 2) intermediate variables are: TVPS, BPNL, COHB
# 3) response (resp) is defined by:

resp <- function(data){
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- resp(smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- cbind(smoking, response)

formula <- TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

mypredict.inclass <- function(object, newdata){
  res <- predict.inclass(object = object, cFUN = resp, newdata = newdata)
  return(res)
}

# Estimation per leave-one-out estimate for the misclassification is 
# 36.36\% (Hand et al., 2001), using indirect classification with 
# linear models

errorest(formula, data = smoking, model = inclass, predict = mypredict.inclass,
         estimator = "cv", iclass = "response", pFUN = lm,
         est.para=list(k=nrow(smoking)))

}

\keyword{misc}
