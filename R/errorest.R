# $Id: errorest.R,v 1.7 2002/03/26 16:29:15 hothorn Exp $

errorest <- function(formula, data, subset, na.action,
                     model=NULL, predict=NULL, iclass=NULL,
                     estimator = c("cv", "boot", "632plus"),
                     est.para = list(k= 10, nboot = 25), ...) {

  m <- match.call(expand.dots = FALSE)
  estimator <- match.arg(estimator)

  # Check if model is "inclass"

  if(paste(m$model) == "inclass") {
    INCLASS <- TRUE
    if (is.null(iclass)) 
      stop("no class membership variable for indirect classification given")
    else
      diagnosis <- iclass
  } else {
    INCLASS <- FALSE
    diagnosis <- attr(terms(formula[-3]), "term.labels")
    if (length(diagnosis) != 1)
      stop("multiple responses not allowed")
  } 
    
  # Check for correct arguments

  if(!INCLASS && class(formula) == "flist") 
    stop("class(formula) = flist is only specified for model = inclass")
  if(INCLASS && is.null(predict)) 
    stop("for model = inclass a prediction model has to be specified")
  if(is.null(model)) 
    stop("no classifier specified")

  # we do not evaluated the formula here but inside "model"

  if (missing(data) & !INCLASS) {
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    m$model <- NULL
    m$predict <- NULL
    m$estimator <- NULL
    m$est.para <- NULL
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    y <- mf[[response]]
    X <- as.data.frame(mf[-response])
    data <- cbind(y, X)
  } else {
    if (INCLASS & missing(data))
      stop("need data argument for indirect classification")
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if(!missing(subset)) data <- as.data.frame(data[subset, ])
  N <- nrow(data)
  y <- data[,diagnosis]
  names(y) <- row.names(data)
  if (!INCLASS)
    DNAME <- paste(diagnosis, "on", attr(terms(formula[-2]), "term.labels"))
  else 
    DNAME <- diagnosis

  # misclassification or mean squared error?

  CLASS <- is.factor(y)

  # k-fold cross-validation

  if(estimator == "cv") {
    if(is.null(est.para$k)) stop("k for k-fold cross-validation is missing")
    k <- est.para$k
    a <- kfoldcv(k, N)
    myindx <- sample(1:N, N)

    cverr <- c()
    for(i in 1:k) {
      if (i > 1)    
        tindx <- myindx[(sum(a[1:(i-1)])+1):sum(a[1:i])]
      else
        tindx <- myindx[1:a[1]]
      X <- as.data.frame(data[-tindx, ])
      if (!missing(na.action))
        mymodel <- model(formula, data = X, na.action = na.action, ...)
      else
        mymodel <- model(formula, data = X, ...)
      if(CLASS)
        pred <- factor(predict(mymodel, newdata = as.data.frame(data[tindx, ])), levels = levels(y))
      if(!CLASS)
        pred <- predict(mymodel, newdata = as.data.frame(data[tindx,]))
      if (CLASS & !is.factor(pred)) 
        stop("predict does not return predicted classes")
      if (length(pred) != length(tindx)) 
        stop("different length of data and prediction")
      if(CLASS) {
        err <- 1 - mean(pred == y[tindx], na.rm = TRUE)
     } else {
        err <- sqrt(mean((pred - y[tindx])^2, na.rm = TRUE))
      }
      cverr <- c(cverr,err)
    }
    err <- mean(cverr)
  }


  if(estimator == "boot" || estimator == "632plus") {
    if(is.null(est.para$nboot)) 
      stop("nboot has to be specified") 
    else 
      nboot <- est.para$nboot
    
    if(CLASS)
      fun <- function(x)
        ifelse(all(is.na(x)), 0, mean(as.integer(x), na.rm = TRUE))
    else
      fun <- function(x)
        ifelse(all(is.na(x)), 0, mean(x, na.rm = TRUE))

    bootindx <- c()

    for(i in 1:nboot) {
      tindx <- sample(1:N, N, replace = TRUE)
      X <- as.data.frame(data[tindx, ])
      if (!missing(na.action))
        mymodel <- model(formula, data = X, na.action = na.action, ...)
      else
        mymodel <- model(formula, data = X,  ...)
      if(CLASS)
        pred <- factor(predict(mymodel, newdata = as.data.frame(data[-tindx, ])), levels = levels(y))
      if(!CLASS)
        pred <- predict(mymodel, newdata = as.data.frame(data[-tindx,]))                    
      if (CLASS & !is.factor(pred)) 
        stop("predict does not return predicted classes")
      if (length(pred) != length(y[-tindx])) 
        stop("different length of data and prediction")

      bootindx <- cbind(bootindx, rep(NA, N))

      if(CLASS)
        bootindx[-tindx, ncol(bootindx)] <- (pred != y[-tindx])
      else
        bootindx[-tindx, ncol(bootindx)] <- sqrt((pred - y[-tindx])^2)
    }
    # mean ( mean deviation (MSE or misclassification rate) for each 
    # observation over the bootstrap samples )
    one <- mean(apply(bootindx, 1, fun))			 
  								 
    err <- one
    expb <- rep(0, nboot)
    for(i in 1:nboot)
      expb[i] <- mean(apply(bootindx[,-i], 1, fun))
    sdint <- sqrt( ((nboot - 1)/nboot)*sum((expb - mean(expb))^2) )
  }

  if(estimator == "632plus") {
    if (!missing(na.action))
      full.model <- model(formula, data = data, na.action = na.action, ...)
    else
      full.model <- model(formula, data = data, ...)
    full.pred <- predict(full.model, newdata = data)
    if (CLASS & !is.factor(full.pred)) 
      stop("predict does not return predicted classes")
    if (length(full.pred) != length(y)) 
      stop("different length of data and prediction")

    if(CLASS)
      resubst <- mean(full.pred != y, na.rm = TRUE)
    else
      stop("632plus implemented for classification problems only")
    
    fun <- function(x, y) ifelse(x==y, 0, 1)
    y <- y[!is.na(y) & !is.na(full.pred)]
    full.pred <- full.pred[!is.na(y) & !is.na(full.pred)]
    gamma <- sum(outer(y, full.pred, fun))/(length(y)^2)

    r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    weight <- .632/(1-.368*r)
    err <- (1-weight)*resubst + weight*one 
    # res <- list(err = err, sdint = sdint) <- gilt hier nicht
  }
  res <- list(err = err, estimator=estimator,
              para=est.para, data.name = DNAME, class = CLASS)
  if(estimator == "boot")
    res <- list(err = err, sd =sdint, estimator=estimator,
                para=est.para, data.name = DNAME, class = CLASS)
  class(res) <- "errorest"
  return(res)
}

print.errorest <- function(x, digits=4, ...)
{
    cat("\n")
    if (x$estimator == "cv")
      method <- paste(x$para$k, "-fold cross-validation", sep="")
    if (x$estimator == "boot")
      method <- "Bootstrap"
    if (x$estimator == "632plus")
      method <- ".632+ Bootstrap"
    method <- paste(method, "estimator of")
    if (x$class)
      method <- paste(method, "misclassification error\n")
    else
      method <- paste(method, "mean squared error\n")
    if (x$estimator == "cv")
      method <- paste(method, "with", x$para$k, "replications")
    else 
      method <- paste(method, "with", x$para$nboot, "bootstrap replications")
  
    cat(method, "\n")
    cat("\n")
    cat("Data: ", x$data.name, "\n")
    cat("\n")
    error <- "Error"
    if (!is.na(x$err))
        cat(error, round(x$err, digits), "\n")
    if (!is.null(x$sdint))
        cat("estimatated standard error", round(x$sdint, digits), "\n")

}
