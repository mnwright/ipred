#$Id: cv.R,v 1.9 2002/12/03 11:41:56 hothorn Exp $

cv <- function(y, ...) UseMethod("cv")

cv.default <- function(y, ...) {
  if (is.numeric(y)) {
    class(y) <- "numeric"
    return(cv(y, ...))
  } else {
    stop(paste("Do not know how to handle objects of class", class(y)))
  }
}

cv.factor <- function(y, X, model, predict, k=10, random=TRUE, 
                      strat=FALSE, predictions=NULL, iformula=NULL, ...) {

  # k-fold cross-validation of misclassification error

  N <- length(y)
  classes <- levels(y)

  if (!is.data.frame(X)) X <- as.data.frame(X)
  if (is.null(k)) k <- 10
  if (is.null(random)) random <- TRUE
  if (is.null(strat)) strat <- FALSE
  if (is.null(predictions)) predictions <- FALSE

  # handle indirect classification with formula w
  # the y,X framework will not work without it
  if (!is.null(iformula)) {
    if (inherits(iformula, "formula") || inherits(iformula, "flist"))
      INCLASS <- TRUE
    else
      stop("iformula must be of class formula or flist") 
  } else { 
    INCLASS <- FALSE
  }

  # to reproduce results, either use `set.seed' or a fixed partition of 
  # the samples
  if (random) 
    myindx <- sample(1:N, N)
  else 
    myindx <- 1:N

  y <- y[myindx]
  X <- X[myindx,]

  # determine an appropriate splitting for the sample size into
  # k roughly equally sized parts

  mysplit <- ssubset(y, k, strat=strat)

  allpred <- vector(mode="character", length=N)
  fu <- function(x) levels(x)[as.integer(x)]
  mydata <- cbind(y,X)
  for(i in 1:k) {
    tindx <- mysplit[[i]]
    if (INCLASS) 
      mymodel <- model(iformula, data=X[-tindx,], ...)
    else
      mymodel <- model(y ~ ., data=mydata[-tindx,], ...)

    # we assume predict to return factor levels
    pred <- factor(predict(mymodel, newdata = X), levels=classes)
    
    # <FIXME>
    # there is no c() for factors which preserves the levels, isn't it?
    # use characters
    allpred[tindx] <- fu(pred[tindx])
    # </FIXME>
  }
  allpred <- factor(allpred, levels=classes)
  allpred <- allpred[order(myindx)]
  err <- mean(allpred != y[order(myindx)], na.rm = TRUE)
  if (predictions)
    RET <- list(error = err, k = k, predictions=allpred)
  else 
    RET <- list(error = err, k = k)
  class(RET) <- "cvclass"
  RET
}

cv.numeric <- function(y, X, model, predict, k=10, random=TRUE,
                       predictions=NULL, strat=NULL, ...) {

  # k-fold cross-validation of mean squared error 

  N <- length(y)

  if (!is.data.frame(X)) X <- as.data.frame(X)
  if (is.null(k)) k <- 10
  if (is.null(random)) random <- TRUE
  if (is.null(predictions)) predictions <- FALSE

  # determine an appropriate splitting for the sample size into
  # k roughly equally sized parts

  a <- kfoldcv(k, N)

  # to reproduce results, either use `set.seed' or a fixed partition of
  # the samples
  if (random)
    myindx <- sample(1:N, N)
  else
    myindx <- 1:N

  allpred <- rep(0, N)
  for(i in 1:k) {
    if (i > 1)
      tindx <- myindx[(sum(a[1:(i-1)])+1):sum(a[1:i])]
    else
      tindx <- myindx[1:a[1]]
    
    mymodel <- model(y ~ ., data=cbind(y, X)[-tindx,], ...)

    pred <- predict(mymodel, newdata = X[tindx,])
    if (!is.numeric(pred)) stop("predict does not return numerical values")
    allpred[tindx] <- pred
  }
  err <- sqrt(mean((pred - y[tindx])^2, na.rm = TRUE))
  allpred <- allpred[order(myindx)]
  if (predictions)
    RET <- list(error = err, k = k, predictions=allpred)
  else
    RET <- list(error = err, k = k)
  class(RET) <- "cvreg" 
  RET  
}

cv.Surv <- function(y, X=NULL, model, predict, k=10, random=TRUE,
                    predictions=FALSE, strat=FALSE, ...) {

  # k-fold cross-validation of Brier's score

  if (is.null(predictions)) predictions <- FALSE
  if(is.null(random)) random <- TRUE
  if (is.null(predictions)) predictions <- FALSE
  if (is.null(strat)) strat <- FALSE


  N <- length(y[,1])
  if(is.null(random)) random <- TRUE
  if(is.null(k)) k <- 10
  if (is.null(X)) X <- rep(1, N)
  # try to find out if model is "survfit"
  options(show.error.messages = FALSE)
  a <- try(model(y))
  if (inherits(a, "survfit")) {
    KM <- TRUE
    if (length(X) != N) stop("only one covariable for survfit allowed")
  } else {
    KM <- FALSE
  }
  options(show.error.messages = TRUE)

  if (!KM)
    if (!is.data.frame(X)) X <- as.data.frame(X)
  
  if(is.null(k)) stop("k for k-fold cross-validation is missing")

  # determine an appropriate splitting for the sample size into
  # k roughly equally sized parts

  a <- kfoldcv(k, N)

  # to reproduce results, either use `set.seed' or a fixed partition of
  # the samples
  if (random)
    myindx <- sample(1:N, N)
  else
    myindx <- 1:N

  cverr <- c()
  for(i in 1:k) {
    if (i > 1)
      tindx <- myindx[(sum(a[1:(i-1)])+1):sum(a[1:i])]
    else
      tindx <- myindx[1:a[1]]

    if (KM)
      mymodel <- survfit(y ~ X, subset=(-tindx), ...)
    else 
      mymodel <- model(y ~ ., data=cbind(y, X)[-tindx,], ...)

    pred <- predict(mymodel, newdata = as.data.frame(X[tindx, ]))
    if (is.list(pred)) {
      if (!inherits(pred[[1]], "survfit") && !inherits(pred, "survfit"))
        stop("predict does not return a list of survfit objects")
    } else {
      stop("predict does not return a list of survfit objects")
    }

    err <- sbrier(y[tindx], pred)
    cverr <- c(cverr,rep(err, length(tindx)))
  }
  RET <- list(error = mean(cverr), k=k)
  class(RET) <- "cvsurv" 
  RET  
}
       