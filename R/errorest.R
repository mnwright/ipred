# $Id: errorest.R,v 1.19 2003/02/25 15:41:40 hothorn Exp $

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
    # just extract the data.frame, NA handling here
    # make sure to leave the time and censoring variable here
    # for "Surv(time, cens) ~ ." formulas
    # delete terms attribute 
    attr(mf, "terms") <- NULL
    y <- mf[,response]
    if (!NOPRED & !is.Surv(y))
      data <- mf
    else
      data <- data[complete.cases(data),]

    estimator <- match.arg(estimator)

    if(is.null(model)) 
      stop("no model specified")

    switch(estimator, "cv" = {
      RET <- cv(y, formula, data, model=model, predict=predict, 
                k=est.para$k, random=est.para$random,
                predictions=est.para$predictions, strat=est.para$strat, ...)
    }, "boot" = {
      RET <- bootest(y, formula, data, model=model, predict=predict,
                     nboot=est.para$nboot, ...)
    }, "632plus" = {
      RET <- bootest(y, formula, data, model=model, predict=predict,
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
    RET <- cv(y, formula, data=X, model=model, predict=predict,
              k=est.para$k, random=est.para$random, ...)
    }, "boot" = {
      RET <- bootest(y, formula, data=X, model=model, predict=predict,
                     nboot=est.para$nboot, ...)
    }, "632plus" = {
      RET <- bootest(y, formula, data=X, model=model, predict=predict,
                     nboot=est.para$nboot, bc632plus=TRUE, ...)
  })
  RET
}
