\name{rsurv}
\alias{rsurv}
\title{ Simulate Survival Data }
\description{
  Simulation Setup for Survival Data.
}
\usage{
rsurv(N, model=c("A", "B", "C", "D", "tree"), gamma=NULL, fact=1, pnon=10,
      gethaz=FALSE)
}
\arguments{
  \item{N}{ number of observations. }
  \item{model}{ type of model. }
  \item{gamma}{simulate censoring time as runif(N, 0, gamma). Defaults to
               \code{NULL} (no censoring).}
  \item{fact}{scale parameter for \code{model=tree}.}
  \item{pnon}{number of additional non-informative variables for the tree
              model.}
  \item{gethaz}{logical, indicating wheather the hazard rate for each 
                observation should be returned.}
}
\details{
  Simulation setup similar to configurations used in LeBlanc and Crowley
(1992) or Keles and Segal (2002) as well as a tree model used in Hothorn et
al. (2002). See Hothorn et al. (2003) for the details. 
}
\value{
  A data  frame with elements \code{time}, \code{cens}, \code{X1} ...
\code{X5}. If \code{pnon} > 0, additional noninformative covariables are
added. If \code{gethaz=TRUE}, the \code{hazard} attribute returns the hazard
rates.
}
\references{ 

 M. LeBlanc and J. Crowley (1992), Relative Risk Trees for 
 Censored Survival Data. \emph{Biometrics} \bold{48}, 411--425.

 S. Keles and M. R. Segal (2002), Residual-based tree-structured 
 survival analysis. \emph{Statistics in Medicine}, \bold{21}, 313--326.

 Torsten Hothorn, Berthold Lausen, Axel Benner and Martin
 Radespiel-Troeger (2004), Bagging Survival Trees.
 \emph{Statistics in Medicine}, \bold{23}(1), 77--91.
 

}
\examples{

# 3*X1 + X2
simdat <- rsurv(500, model="C")
coxph(Surv(time, cens) ~ ., data=simdat)

}
\keyword{survival}
