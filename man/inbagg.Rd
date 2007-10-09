\name{inbagg}
\alias{inbagg}
\alias{inbagg.default}
\alias{inbagg.data.frame}

\title{Indirect Bagging}

\description{
 Function to perform the indirect bagging and subagging.
}

\usage{
\method{inbagg}{data.frame}(formula, data, pFUN=NULL, 
  cFUN=list(model = NULL, predict = NULL, training.set = NULL), 
  nbagg = 25, ns = 0.5, replace = FALSE, ...)
}

\arguments{
  \item{formula}{formula. A \code{formula} specified as \code{y~w1+w2+w3~x1+x2+x3} describes how to model the intermediate variables \code{w1, w2, w3} and the response variable \code{y}, if no other formula is specified by the elements of \code{pFUN} or in \code{cFUN}}
  \item{data}{data frame of explanatory, intermediate and response variables.}
  \item{pFUN}{list of lists, which describe models for the intermediate variables, details are given below.}
  \item{cFUN}{either a fixed function with argument \code{newdata} and returning the class membership by default, or a list specifying a classifying model, similar to one element of \code{pFUN}. Details are given below.}
  \item{nbagg}{number of bootstrap samples.}
  \item{ns}{proportion of sample to be drawn from the learning
    sample. By default, subagging with 50\% is performed, i.e. draw
    0.5*n out of n without replacement.}
  \item{replace}{logical. Draw with or without replacement.}
  \item{\dots}{additional arguments (e.g. \code{subset}).}
}
\details{
 A given data set is subdivided into three types of variables: explanatory, intermediate and response variables.\cr

 Here, each specified intermediate variable is modelled separately
 following \code{pFUN}, a list of lists with elements specifying an
arbitrary number of models for the intermediate variables and an
optional element \code{training.set = c("oob", "bag", "all")}. The
element \code{training.set} determines whether, predictive models for
the intermediate are calculated based on the out-of-bag sample
(\code{"oob"}), the default, on the bag sample (\code{"bag"}) or on all
available observations (\code{"all"}). The elements of \code{pFUN},
specifying the models for the intermediate variables are lists as
described in \code{\link{inclass}}.
Note that, if no formula is given in these elements, the functional
relationship of \code{formula} is used.\cr

The response variable is modelled following \code{cFUN}.
This can either be a fixed classifying function as described in Peters
et al. (2003) or a list,
which specifies the  modelling technique to be applied. The list
contains the arguments \code{model} (which model to be fitted),
\code{predict} (optional, how to predict), \code{formula} (optional, of
type \code{y~w1+w2+w3+x1+x2} determines the variables the classifying
function is based on) and the optional argument \code{training.set =
  c("fitted.bag", "original", "fitted.subset")}
specifying whether the classifying function is trained on the predicted
observations of the bag sample (\code{"fitted.bag"}),
on the original observations (\code{"original"}) or on the
predicted observations not included in a defined subset
(\code{"fitted.subset"}). Per default the formula specified in
\code{formula} determines the variables, the classifying function is
based on.\cr

Note that the default of \code{cFUN = list(model = NULL, training.set = "fitted.bag")}
uses the function \code{\link[rpart]{rpart}} and
the predict function \code{predict(object, newdata, type = "class")}.
}
\value{
  An object of class \code{"inbagg"}, that is a list with elements
  \item{mtrees}{a list of length \code{nbagg}, describing the prediction
    models corresponding
    to each bootstrap sample. Each element of \code{mtrees}
    is a list with elements \code{bindx} (observations of bag sample),
    \code{btree} (classifying function of bag sample) and \code{bfct} (predictive models for intermediates of bag sample).}
  \item{y}{vector of response values.}
  \item{W}{data frame of intermediate variables.}
  \item{X}{data frame of explanatory variables.}
}

\references{
David J. Hand, Hua Gui Li, Niall M. Adams (2001),
Supervised classification with structured class definitions.
\emph{Computational Statistics \& Data Analysis} \bold{36},
209--225.

Andrea Peters, Berthold Lausen, Georg Michelson and Olaf Gefeller (2003),
Diagnosis of glaucoma by indirect classifiers.
\emph{Methods of Information in Medicine} \bold{1}, 99-103.
}

\author{Andrea Peters <Peters.Andrea@imbe.imed.uni-erlangen.de>}

\seealso{\code{\link[rpart]{rpart}}, \code{\link{bagging}},
\code{\link{lm}}}


\examples{
y <- as.factor(sample(1:2, 100, replace = TRUE))
W <- mvrnorm(n = 200, mu = rep(0, 3), Sigma = diag(3))
X <- mvrnorm(n = 200, mu = rep(2, 3), Sigma = diag(3))
colnames(W) <- c("w1", "w2", "w3") 
colnames(X) <- c("x1", "x2", "x3") 
DATA <- data.frame(y, W, X)


pFUN <- list(list(formula = w1~x1+x2, model = lm, predict = mypredict.lm),
list(model = rpart))

inbagg(y~w1+w2+w3~x1+x2+x3, data = DATA, pFUN = pFUN)
}
\keyword{misc}
