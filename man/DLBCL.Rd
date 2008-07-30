\name{DLBCL}
\alias{DLBCL}
\non_function{}
\title{ Diffuse Large B-Cell Lymphoma }
\usage{data(DLBCL)}
\description{
A data frame with gene expression data from diffuse large B-cell
lymphoma (DLBCL) patients.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{DLCL.Sample}{DLBCL identifier.}
    \item{Gene.Expression}{Gene expression group.}
    \item{time}{survival time in month.}
    \item{cens}{censoring: 0 censored, 1 dead.}
    \item{IPI}{International prognostic index.}
    \item{MGEc.1}{mean gene expression in cluster 1.}
    \item{MGEc.2}{mean gene expression in cluster 2.}
    \item{MGEc.3}{mean gene expression in cluster 3.}
    \item{MGEc.4}{mean gene expression in cluster 4.}
    \item{MGEc.5}{mean gene expression in cluster 5.}
    \item{MGEc.6}{mean gene expression in cluster 6.}
    \item{MGEc.7}{mean gene expression in cluster 7.}
    \item{MGEc.8}{mean gene expression in cluster 8.}
    \item{MGEc.9}{mean gene expression in cluster 9.}
    \item{MGEc.10}{mean gene expression in cluster 10.}
  }
}
\source{
Except of \code{MGE}, the data is published at
\url{http://llmpp.nih.gov/lymphoma/data.shtml}. \code{MGEc.*} is the mean of
the gene expression in each of ten clusters derived by agglomerative average
linkage hierarchical cluster analysis (Hothorn et al., 2002).

}
\references{
Ash A. Alizadeh et. al (2000), Distinct types of diffuse large
B-cell lymphoma identified by gene
expression profiling. \emph{Nature}, \bold{403}, 504--509.

Torsten Hothorn, Berthold Lausen, Axel Benner and Martin
Radespiel-Troeger (2004), Bagging Survival Trees. 
\emph{Statistics in Medicine}, \bold{23}, 77--91.


}
\examples{
data(DLBCL)
survfit(Surv(time, cens), data=DLBCL)

}
\keyword{datasets}
