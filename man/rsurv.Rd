\name{rsurv}
\alias{rsurv}
\title{ Simulate Survival Data }
\description{
  Simulation Setup for Survival Data.
}
\usage{
rsurv(N, model=c("A", "B", "C", "tree"), gamma=1, fact=1, pnon=10)
}
\arguments{
  \item{N}{ number of observations. }
  \item{model}{ type of model. }
  \item{gamma}{simulate censoring time as runif(N, 0, gamma).}
  \item{fact}{scale parameter for \code{model=tree}.}
  \item{pnon}{number of additional non-informative variables for the tree
              model.}
}
\details{
  Simulation setup similar to configurations used in LeBlanc and Crowley
(1992) or Keles and Segal (2002) as well as a tree model used in Hothorn et
al. (2002). See Hothorn et al. (2002) for the details. 
}
\value{
  A data  frame with elements \code{time}, \code{cens}, \code{X1} ...
\code{X5}. If \code{pnon} > 0, additional noninformative covariables are
added.
}
\references{ 

 M. LeBlanc and J. Crowley (1992), Relative Risk Trees for 
 Censored Survival Data. \emph{Biometrics} \bold{48}, 411--425.

 S. Keles and M. R. Segal (2002), Residual-based tree-structured 
 survival analysis. \emph{Statistics in Medicine}, \bold{21}, 313--326.

 T. Hothorn, B. Lausen, A. Benner and M. Radespiel-Troeger (2002), 
 Bagging Survival Trees. \emph{submitted}. Preprint available from
 \url{http://www.mathpreprints.com/math/Preprint/blausen/20020518/2}.

}
\examples{

# 3*X1 + X2
simdat <- rsurv(500, model="C")
coxph(Surv(time, cens) ~ ., data=simdat)

}
\keyword{survival}
