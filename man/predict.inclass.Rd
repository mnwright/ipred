\name{predict.inclass}
\alias{predict.inclass}


\title{Predictions from an Inclass Object}

\description{
Predicts the class membership of new observations through indirect classification.
}

\usage{
\method{predict}{inclass}(object, cFUN, intbag=NULL, newdata, ...)
}

\arguments{
  \item{object}{ object of class \code{inclass}, see \code{\link{inclass}}.}
  \item{cFUN}{ fixed function which classifies the predicted intermediates.}
  \item{intbag}{logical. If \code{pFUN = bagging} bagging either stabilises
the prediction of the intermediate variables or the classifications. If
\code{intbag = TRUE} intermediates are aggregated, bagging of the whole 
classifying procedure is performed otherwise.}
  \item{newdata}{data frame to be classified.}
  \item{\dots}{additional arguments corresponding to the predictive models 
specified in \code{\link{inclass}}. If the intermediate is nominal and
\code{pFUN = rpart} or \code{bagging} one argument has to be \code{type = "class"}. }
}

\details{
Predictions of class memberships are calculated. i.e. values of the intermediate variables are predicted and classified following \code{cFUN}, see \code{\link{inclass}}.
}

\value{
The vector of predicted classes is returned. 
}

\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics \& Data Analysis} \bold{36},
209--225.
}

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\seealso{\code{\link{inclass}}}

\examples{
# Simulation model, classification rule following Hand et al. (2001)

theta90 <- varset(N = 1000, sigma = 0.1, theta = 90, threshold = 0)

dataset <- as.data.frame(cbind(theta90$explanatory, theta90$intermediate))
names(dataset) <- c(colnames(theta90$explanatory), colnames(theta90$intermediate))

classify <- function(Y, threshold = 0) {
  z <- (Y > threshold)
  resp <- as.factor(ifelse((z[,1] + z[,2]) > 1, 1, 0))
  return(resp)
}

formula <- flist(y1+y2~x1+x2)

fit <- inclass(formula, pFUN = lm, data = dataset)

predict.inclass(object = fit, cFUN = classify, newdata = dataset)


data(Smoking)

# explanatory variables are: TarY, NicY, COY, Sex, Age
# intermediate variables are: TVPS, BPNL, COHB
# reponse is defined by:

classify <- function(data){
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- classify(Smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- cbind(Smoking, response)

formula <- TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

fit <- inclass(formula, pFUN = lm, data = smoking)

predict.inclass(object = fit, cFUN = classify, newdata = smoking)


data(GlaucomaMVF)
glaucoma <- GlaucomaMVF[,-67]
# explanatory variables are derived by laser scanning image and intra occular pressure
# intermediate variables are: clv, cs, lora
# response is defined by

classify <- function (data) {
  attach(data) 
  res <- ifelse((!is.na(clv) & !is.na(lora) & clv >= 5.1 & lora >= 
        49.23372) | (!is.na(clv) & !is.na(lora) & !is.na(cs) & 
        clv < 5.1 & lora >= 58.55409 & cs < 1.405) | (is.na(clv) & 
        !is.na(lora) & !is.na(cs) & lora >= 58.55409 & cs < 1.405) | 
        (!is.na(clv) & is.na(lora) & cs < 1.405), 0, 1)
  detach(data)
  factor (res, labels = c("glaucoma", "normal"))
}

fit <- inclass(clv+lora+cs~., data = glaucoma, pFUN = bagging)
predict.inclass(object = fit, cFUN = classify, newdata = glaucoma)

}

\keyword{misc}


