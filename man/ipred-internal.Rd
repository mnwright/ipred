\name{ipred-internal}
\alias{uwhich.max}
\alias{predict.sclass}
\alias{predict.sreg}
\alias{getpartition}
\alias{workhorse}
\alias{getsurv}
\alias{ssubset}
\alias{irpart}
\alias{predict.irpart}
\alias{checkfunArgs}
\title{Internal ipred functions}
\description{
 Internal ipred functions.
}
\usage{
uwhich.max(x)
predict.sclass(object, newdata=NULL, type=c("class", "prob"), \dots)
predict.sreg(object, newdata=NULL, \dots)
getpartition(object, newdata=NULL)
workhorse(y, X, control, comb, bcontrol, thisclass, \dots)
getsurv(obj, times)
ssubset(y, k, strat=TRUE)
irpart(formula, data=NULL, weights, subset,
na.action=na.rpart, method, model=FALSE, x=FALSE, y=TRUE,
parms, control, cost, bcontrol, ...)
predict.irpart(object, newdata = list(),
type = c("vector", "prob", "class", "matrix"), ...)
checkfunArgs(fun, type=c("model", "predict"))
}
\details{
  This functions are not to be called by the user.
}
\keyword{internal}
