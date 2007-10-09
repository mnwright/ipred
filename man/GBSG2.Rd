\name{GBSG2}
\alias{GBSG2}
\docType{data} 
\title{ German Breast Cancer Study Group 2 }
\usage{data(GBSG2)}
\description{
  A data frame containing the observations from the  GBSG2 study.
}
\format{
  This data frame contains the observations of 686 women:
  \describe{
    \item{horTh}{hormonal therapy, a factor at two levels \code{no} and
\code{yes}.}
    \item{age}{of the patients in years.}
    \item{menostat}{menopausal status, a factor at two levels \code{pre} 
(premenopausal) and \code{post} (postmenopausal).}
    \item{tsize}{tumor size (in mm).}
    \item{tgrade}{tumor grade, a ordered factor at levels \code{I < II <
III}.}
    \item{pnodes}{number of positive nodes.}
    \item{progrec}{progesterone receptor (in fmol).}
    \item{estrec}{estrogen receptor (in fmol).}
    \item{time}{recurrence free survival time (in days).}
    \item{cens}{censoring indicator (0- censored, 1- event).}
  }
}
\source{
  \url{http://www.blackwellpublishers.com/rss/Volumes/A162p1.htm}
}
\references{
  M. Schumacher, G. Basert, H. Bojar,  K. Huebner, M. Olschewski, 
  W. Sauerbrei, C. Schmoor, C. Beyerle, R.L.A. Neumann and H.F. Rauschecker
  for the German Breast Cancer Study Group (1994),
  Randomized \eqn{2\times2} trial evaluating hormonal treatment
  and the duration of chemotherapy in node-positive breast cancer patients.
  \emph{Journal of Clinical Oncology}, \bold{12}, 2086--2093.

  W. Sauerbrei and P. Royston (1999). Building multivariable prognostic
  and diagnostic models: transformation of the predictors by using 
  fractional polynomials. \emph{Journal of the Royal Statistics Society
  Series A}, Volume \bold{162}(1), 71--94.
}
\examples{
data(GBSG2)

thsum <- function(x) {
  ret <- c(median(x), quantile(x, 0.25), quantile(x,0.75))
  names(ret)[1] <- "Median"
  ret
}

t(apply(GBSG2[,c("age", "tsize", "pnodes", 
                 "progrec", "estrec")], 2, thsum))

table(GBSG2$menostat)
table(GBSG2$tgrade)
table(GBSG2$horTh)

# pooled Kaplan-Meier

mod <- survfit(Surv(time, cens), data=GBSG2)
# integrated Brier score
sbrier(Surv(GBSG2$time, GBSG2$cens), mod)
# Brier score at 5 years
sbrier(Surv(GBSG2$time, GBSG2$cens), mod, btime=1825)

# Nottingham prognostic index

GBSG2 <- GBSG2[order(GBSG2$time),]

NPI <- 0.2*GBSG2$tsize/10 + 1 + as.integer(GBSG2$tgrade)
NPI[NPI < 3.4] <- 1
NPI[NPI >= 3.4 & NPI <=5.4] <- 2
NPI[NPI > 5.4] <- 3

mod <- survfit(Surv(time, cens) ~ NPI, data=GBSG2)
plot(mod)

pred <- c()
survs <- c()
for (i in sort(unique(NPI)))
    survs <- c(survs, getsurv(mod[i], 1825))

for (i in 1:nrow(GBSG2))
   pred <- c(pred, survs[NPI[i]])

# Brier score of NPI at t=5 years
sbrier(Surv(GBSG2$time, GBSG2$cens), pred, btime=1825)


}
\keyword{datasets}
