\name{sbrier}
\alias{sbrier}
\title{ Model Fit for Survival Data }
\description{
  Model fit for survival data: the integrated Brier score for censored
observations.
}
\usage{
sbrier(obj, pred, btime= range(obj[,1]))
}
\arguments{
  \item{obj}{an object of class \code{Surv}.}
  \item{pred}{predicted values. Either a probability or a list of
               \code{survfit} objects. }
  \item{btime}{numeric vector of times, the integrated Brier score is
               computed if this is of \code{length > 1}. 
               The Brier score at \code{btime}
               is returned otherwise.}
}
\details{
  There is no obvious criterion of model fit for censored data. The Brier
score for censoring as well as it's integrated version were suggested by
Graf et al (1999).

The integrated Brier score is always computed over a subset of the
interval given by the range of the time slot of the survival object \code{obj}.

}
\value{
  The (integrated) Brier score with attribute \code{time} is returned. 
}
\seealso{
  More measures for the validation of predicted surival probabilities 
  are implemented in package \code{pec}.
}
\references{ 

Erika Graf, Claudia Schmoor, Willi Sauerbrei and Martin Schumacher (1999),
Assessment and comparison of prognostic classification schemes for
survival data. \emph{Statistics in Medicine} \bold{18}(17-18), 2529--2545.

}

\examples{

data(DLBCL)
smod <- Surv(DLBCL$time, DLBCL$cens)

KM <- survfit(smod ~ 1)
# integrated Brier score up to max(DLBCL$time)
sbrier(smod, KM)

# integrated Brier score up to time=50
sbrier(smod, KM, btime=c(0, 50))

# Brier score for time=50
sbrier(smod, KM, btime=50)

# a "real" model: one single survival tree with Intern. Prognostic Index
# and mean gene expression in the first cluster as predictors
mod <- bagging(Surv(time, cens) ~ MGEc.1 + IPI, data=DLBCL, nbagg=1)

# this is a list of survfit objects (==KM-curves), one for each observation
# in DLBCL
pred <- predict(mod, newdata=DLBCL)

# integrated Brier score up to max(time)
sbrier(smod, pred)

# Brier score at time=50
sbrier(smod, pred, btime=50)
# artificial examples and illustrations

cleans <- function(x) { attr(x, "time") <- NULL; names(x) <- NULL; x }

n <- 100
time <- rpois(n, 20)
cens <- rep(1, n)

# checks, Graf et al. page 2536, no censoring at all!
# no information: \pi(t) = 0.5 

a <- sbrier(Surv(time, cens), rep(0.5, n), time[50])
stopifnot(all.equal(cleans(a),0.25))

# some information: \pi(t) = S(t)

n <- 100
time <- 1:100
mod <- survfit(Surv(time, cens) ~ 1)
a <- sbrier(Surv(time, cens), rep(list(mod), n))
mymin <- mod$surv * (1 - mod$surv)
stopifnot(all.equal(cleans(a),sum(mymin)/max(time)))

# independent of ordering
rand <- sample(1:100)
b <- sbrier(Surv(time, cens)[rand], rep(list(mod), n)[rand])
stopifnot(all.equal(cleans(a), cleans(b)))

\testonly{
  # total information: \pi(t | X) known for every obs

  time <- 1:10
  cens <- rep(1,10)
  pred <- diag(10)
  pred[upper.tri(pred)] <- 1
  diag(pred) <- 0
  # <FIXME>
  # a <- sbrier(Surv(time, cens), pred)
  # stopifnot(all.equal(a, 0))
  # </FIXME>
}

# 2 groups at different risk

time <- c(1:10, 21:30)
strata <- c(rep(1, 10), rep(2, 10))
cens <- rep(1, length(time))

# no information about the groups

a <- sbrier(Surv(time, cens), survfit(Surv(time, cens) ~ 1))
b <- sbrier(Surv(time, cens), rep(list(survfit(Surv(time, cens) ~1)), 20))
stopifnot(all.equal(a, b))

# risk groups known

mod <- survfit(Surv(time, cens) ~ strata)
b <- sbrier(Surv(time, cens), c(rep(list(mod[1]), 10), rep(list(mod[2]), 10)))
stopifnot(a > b)

}
\keyword{survival}
