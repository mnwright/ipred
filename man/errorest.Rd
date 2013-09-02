\name{errorest}
\alias{errorest}
\alias{errorest.data.frame}
\alias{errorest.default}
\title{ Estimators of Prediction Error }
\description{
Resampling based estimates of prediction error: misclassification error, 
root mean squared error or Brier score for survival data.
}
\usage{
\method{errorest}{data.frame}(formula, data, subset, na.action=na.omit, 
         model=NULL, predict=NULL,
         estimator=c("cv", "boot", "632plus"), 
         est.para=control.errorest(), ...)
}

\arguments{
  \item{formula}{a formula of the form \code{lhs ~ rhs}. 
                 Either describing the model of explanatory and 
                 response variables in the usual way (see \code{\link{lm}}) 
                 or the model between explanatory and intermediate variables
                 in the framework of indirect classification, 
                 see \code{\link{inclass}}.}
  \item{data}{a data frame containing the variables in the model formula 
              and additionally the class membership variable 
              if \code{model = inclass}. \code{data} is required for
              indirect classification, otherwise \code{formula} is evaluated
              in the calling environment.}
  \item{subset}{optional vector, specifying a subset of observations to 
                be used.}
  \item{na.action}{function which indicates what should happen when the data
                   contains \code{NA}'s, defaults to \code{\link{na.omit}}.} 
  \item{model}{function. Modelling technique whose error rate is to be 
               estimated. The function \code{model} can either return an 
               object representing a fitted model or a function with
               argument \code{newdata} which returns predicted values. In
               this case, the \code{predict} argument to \code{errorest} is
               ignored.}
  \item{predict}{function. Prediction method to be used. The vector of 
                 predicted values must have the same length as the the 
                 number of to-be-predicted observations. Predictions 
                 corresponding to missing data must be replaced by \code{NA}.
                 Additionally, \code{predict} has to return predicted values 
                 comparable to the responses (that is: factors for 
                 classification problems). See the example on how to make
                 this sure for any predictor.}
  \item{estimator}{estimator of the misclassification error: 
                   \code{cv} cross-validation, \code{boot} bootstrap or 
                   \code{632plus} bias corrected bootstrap (classification
                   only). }
  \item{est.para}{a list of additional parameters that control the
                  calculation of the estimator, see 
                  \code{\link{control.errorest}} for details.}
  \item{\dots}{additional parameters to \code{model}.}
}
\details{
  The prediction error for classification and regression models as well as
predictive models for censored data using cross-validation or the 
bootstrap can be computed by \code{errorest}. For classification problems,
the estimated misclassification error is returned. The root mean squared
error is computed for regression problems and the Brier score for censored
data (Graf et al., 1999) is reported if the response is censored. 

Any model can be specified as long as it is a function with arguments
\code{model(formula, data, subset, na.action, ...)}. If 
a method \code{predict.model(object, newdata, ...)} is available,
\code{predict} does not need to be specified. However, \code{predict} 
has to return predicted values in the same order and of the same length
corresponding to the response. See the examples below. 

$k$-fold cross-validation and the usual bootstrap estimator with
\code{est.para$nboot} bootstrap replications can be computed for
all kind of problems. The bias corrected .632+ bootstrap
by Efron and Tibshirani (1997) is available for classification problems
only. Use \code{\link{control.errorest}} to specify additional arguments.

\code{errorest} is a formula based interface to the generic functions 
\code{\link{cv}} or \code{\link{bootest}} which implement methods for
classification, regression and survival problems.
}

\value{
 The class of the object returned depends on the class of the response
variable and the estimator used. In each case, it is a list with an element
\code{error} and additional information. \code{print} methods are available
for the inspection of the results.
}

\references{

Brian D. Ripley (1996), \emph{Pattern Recognition and Neural Networks}.
Cambridge: Cambridge University Press.

Bradley Efron and Robert Tibshirani (1997),
Improvements on Cross-Validation: The .632+ Bootstrap Estimator.
\emph{Journal of the American Statistical Association} \bold{92}(438),
548--560.

Erika Graf, Claudia Schmoor, Willi Sauerbrei and Martin Schumacher (1999), 
Assessment and comparison of prognostic classification schemes for 
survival data. \emph{Statistics in Medicine} \bold{18}(17-18), 2529--2545.

Rosa A. Schiavo and David J. Hand (2000), Ten More Years of Error Rate
Research. \emph{International Statistical Review} \bold{68}(3), 296-310. 

David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised Classification with Structured Class Definitions.
\emph{Computational Statistics & Data Analysis} \bold{36},
209--225.


}
\examples{

# Classification

data("iris")
library("MASS")

# force predict to return class labels only
mypredict.lda <- function(object, newdata)
  predict(object, newdata = newdata)$class

# 10-fold cv of LDA for Iris data
errorest(Species ~ ., data=iris, model=lda, 
         estimator = "cv", predict= mypredict.lda)

data("PimaIndiansDiabetes", package = "mlbench")
\dontrun{
# 632+ bootstrap of LDA for Diabetes data
errorest(diabetes ~ ., data=PimaIndiansDiabetes, model=lda,
         estimator = "632plus", predict= mypredict.lda)
}

#cv of a fixed partition of the data
list.tindx <- list(1:100, 101:200, 201:300, 301:400, 401:500,
        501:600, 601:700, 701:768)

errorest(diabetes ~ ., data=PimaIndiansDiabetes, model=lda,
          estimator = "cv", predict = mypredict.lda,
          est.para = control.errorest(list.tindx = list.tindx))

\dontrun{
#both bootstrap estimations based on fixed partitions

list.tindx <- vector(mode = "list", length = 25)
for(i in 1:25) {
  list.tindx[[i]] <- sample(1:768, 768, TRUE)
}

errorest(diabetes ~ ., data=PimaIndiansDiabetes, model=lda,
          estimator = c("boot", "632plus"), predict= mypredict.lda,
          est.para = control.errorest(list.tindx = list.tindx))

}
data("Glass", package = "mlbench")

# LDA has cross-validated misclassification error of
# 38\% (Ripley, 1996, page 98)

# Pruned trees about 32\% (Ripley, 1996, page 230)

# use stratified sampling here, i.e. preserve the class proportions
errorest(Type ~ ., data=Glass, model=lda, 
         predict=mypredict.lda, est.para=control.errorest(strat=TRUE))

# force predict to return class labels
mypredict.rpart <- function(object, newdata)
  predict(object, newdata = newdata,type="class")

library("rpart")
pruneit <- function(formula, ...)
  prune(rpart(formula, ...), cp =0.01)

errorest(Type ~ ., data=Glass, model=pruneit,
         predict=mypredict.rpart, est.para=control.errorest(strat=TRUE))

# compute sensitivity and specifity for stabilised LDA

data("GlaucomaM", package = "TH.data")

error <- errorest(Class ~ ., data=GlaucomaM, model=slda,
  predict=mypredict.lda, est.para=control.errorest(predictions=TRUE))

# sensitivity 

mean(error$predictions[GlaucomaM$Class == "glaucoma"] == "glaucoma")

# specifity

mean(error$predictions[GlaucomaM$Class == "normal"] == "normal")

# Indirect Classification: Smoking data

data(Smoking)
# Set three groups of variables:
# 1) explanatory variables are: TarY, NicY, COY, Sex, Age
# 2) intermediate variables are: TVPS, BPNL, COHB
# 3) response (resp) is defined by:

resp <- function(data){
  data <- data[, c("TVPS", "BPNL", "COHB")]
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- resp(Smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- cbind(Smoking, response)

formula <- response~TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

# Estimation per leave-one-out estimate for the misclassification is 
# 36.36\% (Hand et al., 2001), using indirect classification with 
# linear models
\dontrun{
errorest(formula, data = smoking, model = inclass,estimator = "cv", 
         pFUN = list(list(model=lm, predict = mypredict.lm)), cFUN = resp,  
         est.para=control.errorest(k=nrow(smoking)))
}

# Regression

data("BostonHousing", package = "mlbench")

# 10-fold cv of lm for Boston Housing data
errorest(medv ~ ., data=BostonHousing, model=lm,
         est.para=control.errorest(random=FALSE))

# the same, with "model" returning a function for prediction
# instead of an object of class "lm"

mylm <- function(formula, data) {
  mod <- lm(formula, data)
  function(newdata) predict(mod, newdata)
}

errorest(medv ~ ., data=BostonHousing, model=mylm,
est.para=control.errorest(random=FALSE))


# Survival data

data("GBSG2", package = "TH.data")
library("survival")

# prediction is fitted Kaplan-Meier
predict.survfit <- function(object, newdata) object

# 5-fold cv of Kaplan-Meier for GBSG2 study
errorest(Surv(time, cens) ~ 1, data=GBSG2, model=survfit,
         predict=predict.survfit, est.para=control.errorest(k=5))


}

\keyword{misc}
