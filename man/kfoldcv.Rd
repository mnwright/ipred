\name{kfoldcv}
\alias{kfoldcv}
\title{ Subsamples for k-fold Cross-Validation }
\description{
  Computes feasible sample sizes for the k groups in k-fold cv if N/k is not
an integer.
}
\usage{
kfoldcv(k, N, nlevel=NULL)
}
\arguments{
  \item{k}{ number of groups. }
  \item{N}{ total sample size. }
  \item{nlevel}{ a vector of sample sizes for stratified sampling.}
}
\details{
  If N/k is not an integer, k-fold cv is not unique. Determine meaningful
  sample sizes.
}
\value{
  A vector of length \code{k}.
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

# 10-fold CV with N = 91

kfoldcv(10, 91)	

\testonly{
k <- sample(5:15, 1)
k
N <- sample(50:150, 1)
N
stopifnot(sum(kfoldcv(k, N)) == N)
}

}
\keyword{misc}
