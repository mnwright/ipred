# $Id: kfoldcv.R,v 1.2 2002/03/02 14:21:15 hothorn Exp $

kfoldcv <- function(k,N) {
  if (k > N) stop("k > N")
  fl <- floor(N/k)
  ce <- ceiling(N/k)
  if (fl == ce) rep(fl, k) else 
  c(rep(ce, round((N/k - fl)*k)), rep(fl, round((1 - (N/k - fl))*k)))
}
