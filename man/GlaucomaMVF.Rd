\name{GlaucomaMVF}
\alias{GlaucomaMVF}
\non_function{}
\title{ Glaucoma Database }
\usage{data(GlaucomaMVF)}
\description{
The \code{GlaucomaMVF} data has 170 observations in two classes. 
66 predictors are derived from a confocal laser scanning image of the 
optic nerve head, from a visual field test, a fundus photography and a 
measurement of the intra occular pressure. 
}
\format{
  This data frame contains the following predictors describing the
  morphology of the optic nerve head, the visual field, the intra 
  occular pressure and a membership variable:
  \describe{
    \item{ag}{area global.}
    \item{at}{area temporal.}
    \item{as}{area superior.}
    \item{an}{area nasal.}
    \item{ai}{area inferior.}
    \item{eag}{effective area global.}
    \item{eat}{effective area temporal.}
    \item{eas}{effective area superior.}
    \item{ean}{effective area nasal.}
    \item{eai}{effective area inferior.}
    \item{abrg}{area below reference global.}
    \item{abrt}{area below reference temporal.}
    \item{abrs}{area below reference superior.}
    \item{abrn}{area below reference nasal.}
    \item{abri}{area below reference inferior.}
    \item{hic}{height in contour.}
    \item{mhcg}{mean height contour global.}
    \item{mhct}{mean height contour temporal.}
    \item{mhcs}{mean height contour superior.}
    \item{mhcn}{mean height contour nasal.}
    \item{mhci}{mean height contour inferior.}
    \item{phcg}{peak height contour.}
    \item{phct}{peak height contour temporal.}
    \item{phcs}{peak height contour superior.}
    \item{phcn}{peak height contour nasal.}
    \item{phci}{peak height contour inferior.}
    \item{hvc}{height variation contour.}
    \item{vbsg}{volume below surface global.}
    \item{vbst}{volume below surface temporal.}
    \item{vbss}{volume below surface superior.}
    \item{vbsn}{volume below surface nasal.}
    \item{vbsi}{volume below surface inferior.}
    \item{vasg}{volume above surface global.}
    \item{vast}{volume above surface temporal.}
    \item{vass}{volume above surface superior.}
    \item{vasn}{volume above surface nasal.}
    \item{vasi}{volume above surface inferior.}
    \item{vbrg}{volume below reference global.}
    \item{vbrt}{volume below reference temporal.}
    \item{vbrs}{volume below reference superior.}
    \item{vbrn}{volume below reference nasal.}
    \item{vbri}{volume below reference inferior.}
    \item{varg}{volume above reference global.}
    \item{vart}{volume above reference temporal.}
    \item{vars}{volume above reference superior.}
    \item{varn}{volume above reference nasal.}
    \item{vari}{volume above reference inferior.}
    \item{mdg}{mean depth global.}
    \item{mdt}{mean depth temporal.}
    \item{mds}{mean depth superior.}
    \item{mdn}{mean depth nasal.}
    \item{mdi}{mean depth inferior.}
    \item{tmg}{third moment global.}
    \item{tmt}{third moment temporal.}
    \item{tms}{third moment superior.}
    \item{tmn}{third moment nasal.}
    \item{tmi}{third moment inferior.}
    \item{mr}{mean radius.}
    \item{rnf}{retinal nerve fiber thickness.}
    \item{mdic}{mean depth in contour.}
    \item{emd}{effective mean depth.}
    \item{mv}{mean variability.}
    \item{tension}{intra occular pressure.}
    \item{clv}{corrected loss variance, variability of the visual field.}
    \item{cs}{contrast sensitivity of the visual field.}
    \item{lora}{loss of rim area, measured by fundus photography.}
    \item{Class}{a factor with levels \code{glaucoma} and \code{normal}.}
  }
}
\details{
Confocal laser images of the eye background are taken with the 
Heidelberg Retina Tomograph and variables 1-62 are derived. 
Most of these variables describe either the area or volume in 
certain parts of the papilla and are measured in 
four sectors (temporal, superior, nasal and inferior) as well 
as for the whole papilla (global). The global measurement is, 
roughly, the sum of the measurements taken in the four sector.

The perimeter `Octopus' measures the visual field variables \code{clv} 
and \code{cs}, stereo optic disks photographs were taken with a 
telecentric fundus camera and \code{lora} is derived.

Observations of both groups are matched by age and sex, 
to prevent for possible confounding. 

}
\source{
  Andrea Peters, Berthold Lausen, Georg Michelson and Olaf Gefeller (2002), 
  Diagnosis of glaucoma by indirect classifiers, 
  \emph{Methods of Information in Medicine}, to appear.
}

\note{
 \code{GLaucomMVF} overlaps in some parts with \code{\link{GlaucomaM}}.
}

\examples{
data(GlaucomaMVF)

response <- function (data) {
  attach(data)
  res <- ifelse((!is.na(clv) & !is.na(lora) & clv >= 5.1 & lora >= 
        49.23372) | (!is.na(clv) & !is.na(lora) & !is.na(cs) & 
        clv < 5.1 & lora >= 58.55409 & cs < 1.405) | (is.na(clv) & 
        !is.na(lora) & !is.na(cs) & lora >= 58.55409 & cs < 1.405) | 
        (!is.na(clv) & is.na(lora) & cs < 1.405), 0, 1)
  detach(data)
  factor (res, labels = c("glaucoma", "normal"))
}

errorest(Class~clv+lora+cs~., data = GlaucomaMVF, model=inclass, 
       estimator="cv", pFUN = list(list(model = rpart)), cFUN = response)

}
\keyword{datasets}
