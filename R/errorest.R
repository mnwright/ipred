# $Id: errorest.R,v 1.14 2002/09/26 15:21:47 peters Exp $

control.errorest <- function(k= 10, nboot = 25, strat=FALSE,
                     random=TRUE, predictions=FALSE) {
  if (k < 1) { 
    warning("k < 1, using k=10")
    k <- 10
  }
  if (nboot < 1) {
    warning("nboot < 1, using nboot=25")
    nboot <- 25
  }
  if (!is.logical(strat)) {
    warning("strat is not a logical, using strat=FALSE")
    strat <- FALSE
  }
  if (!is.logical(random)) {
    warning("random is not a logical, using random=TRUE")
    random <- TRUE
  }
  if (!is.logical(predictions)) {
    warning("predictions is not a logical, using predictions=FALSE")
    predictions <- FALSE
  }

  RET <- list(k=k, nboot=nboot, strat=strat, random=random, 
              predictions=predictions)
  return(RET)
}

errorest <- function(formula, data, ...) UseMethod("errorest", data)

errorest.default <- function(formula, data, ...)
  stop(paste("Do not know how to handle objects of class", class(data)))

errorest.data.frame <- function(formula, data, subset, na.action=na.omit,
                     model=NULL, predict=NULL, iclass=NULL,
                     estimator = c("cv", "boot", "632plus"),
                     est.para = control.errorest(), ...) {

  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if(paste(m$model) == "inclass") {
    RET <- errorestinclass(flist(formula), data=data, subset, na.action,
           model, predict, iclass, estimator, est.para, ...)
    RET$call <- cl
  } else { 
    if(missing(formula)
      || (length(formula) != 3)
      || (length(attr(terms(formula[-3]), "term.labels")) != 1))
    stop("formula missing or incorrect")
    NOPRED <- (length(attr(terms(formula[-2]), "term.labels")) < 1) 
    if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    m$model <- NULL
    m$predict <- NULL
    m$estimator <- NULL
    m$est.para <- NULL

    mf <- eval(m, parent.frame())

    response <- attr(attr(mf, "terms"), "response")
    # just extract the data.frame, no handling of contrasts or NA's here.
    # this is done by rpart or the user supplied methods
    DATA <- list(y = mf[,response], X = mf[,-response]) 
    names(DATA) <- c("y", "X")
    N <- nrow(DATA)

    # y ~ 1 is needed for overall Kaplan-Meier, for example. *.Surv can deal
    # with this
    if (NOPRED)
      DATA$X <- NULL

    estimator <- match.arg(estimator)

    if(is.null(model)) 
      stop("no classifier specified")

    switch(estimator, "cv" = {
      RET <- cv(DATA$y, DATA$X, model=model, predict=predict, 
                k=est.para$k, random=est.para$random,
                predictions=est.para$predictions, strat=est.para$strat, ...)
    }, "boot" = {
      RET <- bootest(DATA$y, DATA$X, model=model, predict=predict,
                     nboot=est.para$nboot, ...)
    }, "632plus" = {
      RET <- bootest(DATA$y, DATA$X, model=model, predict=predict,
                     nboot=est.para$nboot, bc632plus=TRUE, ...)
    })
  }
  RET$call <- cl
  return(RET)
}

errorestinclass <- function(formula, data, subset=NULL, na.action=NULL, 
                     model=NULL, predict=NULL, iclass=NULL,
                     estimator = c("cv", "boot", "632plus"),
                     est.para = control.errorest(), ...) {
  if (is.null(data)) stop("data argument required but not given")
  if (is.null(iclass)) 
    stop("no class membership variable for indirect classification given")
  if (!(iclass %in% colnames(data))) 
    stop("membership variable not in given data")

  # <FIXME> 
#  data <- data[complete.cases(data),]
  # </FIXME>

  iclassindx <- which(colnames(data) == iclass)

  y <- data[,iclassindx]
  if (!is.factor(y)) stop("iclass is not a factor")
  X <- data[,-iclassindx]

  if(is.null(model))
      stop("no classifier specified")

  switch(estimator, "cv" = {
    RET <- cv(y, X, model=model, predict=predict,
              k=est.para$k, random=est.para$random, iformula=formula, ...)
    }, "boot" = {
      RET <- bootest(y, X, model=model, predict=predict,
                     nboot=est.para$nboot, iformula=formula, ...)
    }, "632plus" = {
      RET <- bootest(y, X, model=model, predict=predict,
                     nboot=est.para$nboot, iformula=formula, bc632plus=TRUE, ...)
  })
  RET
}
