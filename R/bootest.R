# $Id: bootest.R,v 1.7 2003/02/04 17:41:13 hothorn Exp $

bootest <- function(y, ...) {
  if(is.null(class(y)))
    class(y) <- data.class(y)
  UseMethod("bootest", y, ...)
}

bootest.default <- function(y, ...) {
  stop(paste("Do not know how to handle objects of class", class(y)))
}



bootest.factor <- function(y, X, model, predict, 
                           nboot=25, bc632plus=FALSE, iformula=NULL, ...) {
  
  # bootstrap estimator of misclassification error

  N <- length(y)
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  classes <- levels(y)
  
  if (!is.data.frame(X)) X <- as.data.frame(X)

  if(is.null(nboot)) stop("number of bootstrap replications is missing")

  # handle indirect classification
  if (!is.null(iformula)) {
    if (inherits(iformula, "formula") || inherits(iformula, "flist"))
      INCLASS <- TRUE 
    else
      stop("iformula must be of class formula or flist")
  } else {
    INCLASS <- FALSE
  }

  mydata <- cbind(y,X)

  for(i in 1:nboot) {
    tindx <- sample(1:N, N, replace = TRUE)
    if (INCLASS)
      mymodel <- model(iformula, data = X[tindx,], ...)
    else
      mymodel <- model(y ~ ., data = mydata[tindx,], ...)

    pred <- factor(predict(mymodel, newdata = X[-tindx, ]), levels=classes)
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")

    bootindx[-tindx, i] <- (pred != y[-tindx])
  }

  fun <- function(x)
        ifelse(all(is.na(x)), 0, mean(as.integer(x), na.rm = TRUE))

  one <- mean(apply(bootindx, 1, fun))

  if (bc632plus) {
    if(INCLASS)
      full.model <- model(iformula, data = X, ...)
    else
      full.model <- model(y ~ ., data = X, ...)
    full.pred <- predict(full.model, newdata = X)
    resubst <- mean(full.pred != y, na.rm = TRUE)

    y <- y[!is.na(y) & !is.na(full.pred)]
    full.pred <- full.pred[!is.na(y) & !is.na(full.pred)]
    gamma <- sum(outer(y, full.pred, function(x, y) ifelse(x==y, 0, 1) ))/
                 (length(y)^2)
    r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    weight <- .632/(1-.368*r)
    err <- (1-weight)*resubst + weight*one
    RET <- list(error = err, nboot=nboot, bc632plus=TRUE)
  } else {
    err <- one
    expb <- rep(0, nboot)
    for(i in 1:nboot)
      expb[i] <- mean(apply(bootindx[,-i], 1, fun))
    sdint <- sqrt( ((nboot - 1)/nboot)*sum((expb - mean(expb))^2) )
    RET <- list(error = err, sd=sdint, bc632plus=FALSE, nboot=nboot)
  }
  class(RET) <- "bootestclass"
  RET
}

bootest.numeric <- function(y, X, model, predict, 
                           nboot=25, bc632plus=FALSE, ...) {
  
  # bootstrap estimator of root of mean squared error 

  if (bc632plus) stop("cannot compute 632+ estimator of mean squared error")

  N <- length(y)
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  
  if (!is.data.frame(X)) X <- as.data.frame(X)

  if(is.null(nboot)) stop("number of bootstrap replications is missing")


  for(i in 1:nboot) {
    tindx <- sample(1:N, N, replace = TRUE)
    mymodel <- model(y ~ ., data = cbind(y, X)[tindx,], ...)
    pred <- predict(mymodel, newdata = X[-tindx, ])
    if (!is.numeric(pred)) stop("predict does not return numerical values")
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")

    bootindx[-tindx, i] <- sqrt((pred - y[-tindx])^2)
  }

  fun <- function(x)
        ifelse(all(is.na(x)), 0, mean(x, na.rm = TRUE))

  err <- mean(apply(bootindx, 1, fun))
  RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestreg"
  RET
}

bootest.Surv <- function(y, X, model, predict, 
                           nboot=25, bc632plus=FALSE, ...) {
  
  # bootstrap estimator of Brier's score

  if (bc632plus) stop("cannot compute 632+ estimator of Brier's score")

  N <- dim(y)[1]
  bootindx <- matrix(NA, ncol=nboot, nrow=N)

  if(is.null(nboot)) stop("number of bootstrap replications is missing")

  if (is.null(X)) X <- rep(1, N)
  # try to find out if model estimates KM
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

  for(i in 1:nboot) {
    tindx <- sample(1:N, N, replace = TRUE)

    if (KM)
      mymodel <- survfit(y ~ X, subset=tindx, ...)
    else 
      mymodel <- model(y ~ ., data=cbind(y, X)[tindx, ], ...)

    pred <- predict(mymodel, newdata = X[-tindx, ])

    if (is.list(pred)) {
      if (!inherits(pred[[1]], "survfit") && !inherits(pred, "survfit"))
        stop("predict does not return a list of survfit objects")
    } else {
      stop("predict does not return a list of survfit objects")
    }

    bootindx[-tindx, i] <- sbrier(y[-tindx], pred)
  }

  fun <- function(x)
        ifelse(all(is.na(x)), 0, mean(x, na.rm = TRUE))

  err <- mean(apply(bootindx, 1, fun))
  RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestsurv"
  RET
}

