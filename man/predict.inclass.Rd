\name{predict.inclass}
\alias{predict.inclass}

\title{Predictions from an Inclass Object}

\description{
Predicts the class membership of new observations through indirect
classification.
}

\usage{
\method{predict}{inclass}(object, newdata, ...)
}

\arguments{
  \item{object}{ object of class \code{inclass}, see \code{\link{inclass}}.}
  \item{newdata}{data frame to be classified.}
  \item{\dots}{additional arguments corresponding to the predictive models 
specified in \code{\link{inclass}}.}
}
\details{
Predictions of class memberships are calculated. i.e. values of the
intermediate variables are predicted and classified following \code{cFUN},
see \code{\link{inclass}}.
}
\value{
The vector of predicted classes is returned.
}
\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics & Data Analysis} \bold{36},
209--225.

Andrea Peters, Berthold Lausen, Georg Michelson and Olaf Gefeller (2003),
Diagnosis of glaucoma by indirect classifiers.
\emph{Methods of Information in Medicine} \bold{1}, 99-103.
}


\seealso{\code{\link{inclass}}}

\examples{
\dontrun{
# Simulation model, classification rule following Hand et al. (2001)

theta90 <- varset(N = 1000, sigma = 0.1, theta = 90, threshold = 0)

dataset <- as.data.frame(cbind(theta90$explanatory, theta90$intermediate))
names(dataset) <- c(colnames(theta90$explanatory),
colnames(theta90$intermediate))

classify <- function(Y, threshold = 0) {
  Y <- Y[,c("y1", "y2")]
  z <- (Y > threshold)
  resp <- as.factor(ifelse((z[,1] + z[,2]) > 1, 1, 0))
  return(resp)
}

formula <- response~y1+y2~x1+x2

fit <- inclass(formula, data = dataset, pFUN = list(list(model = lm)), 
 cFUN = classify)

predict(object = fit, newdata = dataset)


data("Smoking", package = "ipred")

# explanatory variables are: TarY, NicY, COY, Sex, Age
# intermediate variables are: TVPS, BPNL, COHB
# reponse is defined by:

classify <- function(data){
  data <- data[,c("TVPS", "BPNL", "COHB")]
  res <- t(t(data) > c(4438, 232.5, 58))
  res <- as.factor(ifelse(apply(res, 1, sum) > 2, 1, 0))
  res
}

response <- classify(Smoking[ ,c("TVPS", "BPNL", "COHB")])
smoking <- cbind(Smoking, response)

formula <- response~TVPS+BPNL+COHB~TarY+NicY+COY+Sex+Age

fit <- inclass(formula, data = smoking, 
  pFUN = list(list(model = lm)), cFUN = classify)


predict(object = fit, newdata = smoking)
}

data("GlaucomaMVF", package = "ipred")
library("rpart")
glaucoma <- GlaucomaMVF[,(names(GlaucomaMVF) != "tension")]
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

fit <- inclass(Class~clv+lora+cs~., data = glaucoma, 
             pFUN = list(list(model = rpart)), cFUN = classify)

data("GlaucomaM", package = "TH.data")
predict(object = fit, newdata = GlaucomaM)

}

\keyword{misc}
