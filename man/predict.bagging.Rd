\name{predict.classbagg}
\alias{predict.classbagg}
\alias{predict.regbagg}
\alias{predict.survbagg}
\title{ Predictions from Bagging Trees }
\description{
  Predict the outcome of a new observation based on multiple trees.  
}
\usage{
predict.classbagg(object, newdata=NULL, type=c("class", "prob"),
                            aggregation=c("majority", "average", "weighted"), \dots)
predict.regbagg(object, newdata=NULL, aggregation=c("average",
                "weighted"), \dots)
predict.survbagg(object, newdata=NULL,\dots)
}
\arguments{
  \item{object}{object of classes \code{classbagg}, \code{regbagg} or
                \code{survbagg}.}
  \item{newdata}{a data frame of new observations. }
  \item{type}{character string denoting the type of predicted value
             returned for classification trees. Either \code{class} 
             (predicted classes are returned) or \code{prob} 
             (estimated class probabilities are returned).}
  \item{aggregation}{character string specifying how to aggregate, see below.}
  \item{...}{additional arguments, currently not passed to any function.}
}
\details{
  There are (at least) three different ways to aggregate the predictions of
bagging classification trees. Most famous is class majority voting
(\code{aggregation="majority"}) where the most frequent class is returned. The
second way is choosing the class with maximal averaged class probability
(\code{aggregation="average"}). The third method is based on the "aggregated learning
sample", introduced by Hothorn et al. (2003) for survival trees.
The prediction of a new observation is the majority class, mean or
Kaplan-Meier curve of all observations from the learning sample 
identified by the \code{nbagg} leaves containing the new observation.  
For regression trees, only averaged or weighted predictions are possible. 

By default, the out-of-bag estimate is computed if \code{newdata} is NOT
specified. Therefore, the predictions of \code{predict(object)} are "honest"
in some way (this is not possible for combined models via \code{comb} in
\code{\link{bagging}}). 
If you like to compute the predictions for the learning sample
itself, use \code{newdata} to specify your data. 

}
\value{
  The predicted class or estimated class probabilities are returned for
classification trees. The predicted endpoint is returned in regression
problems and the predicted Kaplan-Meier curve is returned for survival
trees. 
}

\author{ Torsten.Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\references{ 

Leo Breiman (1996), Bagging Predictors. \emph{Machine Learning}
\bold{24}(2), 123--140.

Torsten Hothorn, Berthold Lausen, Axel Benner and Martin
Radespiel-Troeger (2004), Bagging Survival Trees.
\emph{Statistics in Medicine}, \bold{23}(1), 77--91.
 
}

\examples{

data(Ionosphere)
Ionosphere$V2 <- NULL # constant within groups

# nbagg = 10 for performance reasons here
mod <- bagging(Class ~ ., data=Ionosphere)

# out-of-bag estimate

mean(predict(mod) != Ionosphere$Class)

# predictions for the first 10 observations

predict(mod, newdata=Ionosphere[1:10,])

predict(mod, newdata=Ionosphere[1:10,], type="prob")

}
\keyword{tree}
