\name{dystrophy}
\alias{dystrophy}
\non_function{}
\title{Detection of muscular dystrophy carriers.}
\usage{data(dystrophy)}
\description{
The \code{dystrophy} data frame has 209 rows and 10 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{OBS}{numeric. Observation number.}
    \item{HospID}{numeric. Hospital ID number.}
    \item{AGE}{numeric, age in years.}
    \item{M}{numeric. Month of examination.}
    \item{Y}{numeric. Year of examination.}
    \item{CK}{numeric. Serum marker creatine kinase.}
    \item{H}{numeric. Serum marker hemopexin.}
    \item{PK}{numeric. Serum marker pyruvate kinase.}
    \item{LD}{numeric. Serum marker lactate dehydroginase.}
    \item{Class}{factor with levels, \code{carrier} and \code{normal}.}
  }
}

\details{
Duchenne Muscular Dystrophy (DMD) is a genetically transmitted disease,
passed from a mother to her children. Affected female offspring usually suffer
no apparent symptoms, male offspring with the disease die at young age.
Although female carriers have no physical symptoms they tend to exhibit
elevated levels of certain serum enzymes or proteins.
\cr
The dystrophy dataset contains 209 observations of 75 female DMD carriers and
134 female DMD non-carrier. It includes 6 variables describing age of the
female and the serum parameters serum marker creatine kinase (CK), serum marker
 hemopexin (H), serum marker pyruvate kinase (PK) and serum marker lactate
dehydroginase (LD). The serum markers CK and H may be measured rather
inexpensive from frozen serum, PK and LD requires fresh serum.
}

\source{
D.Andrews and A. Herzberg (1985), Data. Berlin: Springer-Verlag.
}

\references{
Robert Tibshirani and Geoffry Hinton (1998), Coaching variables for regression and classification. Statistics and Computing 8, 25-33.
}

\examples{
\dontrun{

data("dystrophy")
library("rpart")
errorest(Class~CK+H~AGE+PK+LD, data = dystrophy, model = inbagg, 
pFUN = list(list(model = lm, predict = mypredict.lm), list(model = rpart)), 
ns = 0.75, estimator = "cv")
}
}

\keyword{datasets}
