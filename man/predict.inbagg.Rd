\name{predict.inbagg}
\alias{predict.inbagg}

\title{Predictions from an Inbagg Object}
\description{
Predicts the class membership of new observations through indirect
bagging.
}
\usage{
predict.inbagg(object, newdata, ...)
}
\arguments{
  \item{object}{object of class \code{inbagg}, see \code{\link{inbagg}}.}
  \item{newdata}{data frame to be classified.}
  \item{\dots}{additional argumends corresponding to the predictive models.}
}
\details{
Predictions of class memberships are calculated. i.e. values of the
intermediate variables are predicted following \code{pFUN} and classified following \code{cFUN},
see \code{\link{inbagg}}.
}
\value{
The vector of predicted classes is returned.
}

\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics \& Data Analysis} \bold{36}, 
209--225.

Andrea Peters, Berthold Lausen, Georg Michelson and Olaf Gefeller (2003),
Diagnosis of glaucoma by indirect classifiers.
\emph{Methods of Information in Medicine} \bold{1}, 99-103.
}        

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\seealso{\code{\link{inbagg}}}

\examples{
library(mvtnorm)
y <- as.factor(sample(1:2, 100, replace = TRUE))
W <- mvrnorm(n = 200, mu = rep(0, 3), Sigma = diag(3)) 
X <- mvrnorm(n = 200, mu = rep(2, 3), Sigma = diag(3))
colnames(W) <- c("w1", "w2", "w3")
colnames(X) <- c("x1", "x2", "x3")
DATA <- data.frame(y, W, X)

pFUN <- list(list(formula = w1~x1+x2, model = lm),
list(model = rpart))

RES <- inbagg(y~w1+w2+w3~x1+x2+x3, data = DATA, pFUN = pFUN)
predict(RES, newdata = X)
}
\keyword{misc}
