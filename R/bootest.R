# $Id: bootest.R,v 1.14 2003/03/12 17:06:40 hothorn Exp $

bootest <- function(y, ...) {
  if(is.null(class(y)))
    class(y) <- data.class(y)
  UseMethod("bootest", y, ...)
}

bootest.default <- function(y, ...) {
  stop(paste("Do not know how to handle objects of class", class(y)))
}

bootest.integer <- function(y, ...) {
  bootest.numeric(y, ...)
}

bootest.factor <- function(y, formula, data, model, predict, 
                           nboot=25, bc632plus=FALSE, ...) {
  
  # bootstrap estimator of misclassification error

  N <- length(y)
  nindx <- 1:N
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  classes <- levels(y)
  USEPM <- FALSE
  
  if(!is.data.frame(data)) stop("data is not a data.frame")
  if(nboot <=2) stop("to small number of bootstrap replications")
  if(is.null(nboot)) stop("number of bootstrap replications is missing")

  for(i in 1:nboot) {
    tindx <- sample(nindx, N, replace = TRUE)
    mymodel <- model(formula, data = data[tindx,], ...)

    # check if mymodel is a function which should be used instead of   
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM) 
      pred <- predict(newdata=data)
    else 
      pred <- predict(mymodel, newdata = data)
    if (!is.factor(pred)) stop("predict does not return factor values")
    pred <- factor(pred, levels=classes)[-tindx]
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")

    bootindx[-tindx, i] <- (pred != y[-tindx])
  }

  fun <- function(x)
       ifelse(all(is.na(x)), NA, mean(as.integer(x), na.rm = TRUE))

  one <- mean(apply(bootindx, 1, fun), na.rm = TRUE)

  if (bc632plus) {
    full.model <- model(formula, data = data, ...)
    # check if full.model is a function which should be used instead of
    # predict
    if (is.function(full.model)) {
      predict <- fullmodel
      USEPM <- TRUE
    }

    if (USEPM)
      full.pred <- predict(newdata=data)
    else

    full.pred <- predict(full.model, newdata = data)
    resubst <- mean(full.pred != y, na.rm = TRUE)

    err632 <- 0.368*resubst + 0.632*one

    y <- y[!is.na(y) & !is.na(full.pred)]
    full.pred <- full.pred[!is.na(y) & !is.na(full.pred)]
    gamma <- sum(outer(y, full.pred, function(x, y) ifelse(x==y, 0, 1) ))/
                 (length(y)^2)
    r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    errprime <- min(one, gamma)
    #    weight <- .632/(1-.368*r)
    #    err <- (1-weight)*resubst + weight*one
    err <- err632 + (errprime - resubst)*(0.368*0.632*r)/(1-0.368*r)
    RET <- list(error = err, nboot=nboot, bc632plus=TRUE)
  } else {
    err <- one
    expb <- rep(0, nboot)
    for(i in 1:nboot)
      expb[i] <- mean(apply(bootindx[,-i], 1, fun), na.rm = TRUE)
    sdint <- sqrt( ((nboot - 1)/nboot)*sum((expb - mean(expb))^2) )
    RET <- list(error = err, sd=sdint, bc632plus=FALSE, nboot=nboot)
  }
  class(RET) <- "bootestclass"
  RET
}

bootest.numeric <- function(y, formula, data, model, predict, 
                           nboot=25, bc632plus=FALSE, ...) {
  
  # bootstrap estimator of root of mean squared error 

  if (bc632plus) stop("cannot compute 632+ estimator of mean squared error")
  if (nboot <=2) stop("to small number of bootstrap replications")

  N <- length(y)
  nindx <- 1:N
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  USEPM <- FALSE
  
  if (!is.data.frame(data)) stop("data is not a data.frame")

  if(is.null(nboot)) stop("number of bootstrap replications is missing")


  for(i in 1:nboot) {
    tindx <- sample(nindx, N, replace = TRUE)
    mymodel <- model(formula, data = data[tindx,], ...)
    outbootdata <- subset(data, !(nindx %in% tindx))
    # check if mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM)
      pred <- predict(newdata=outbootdata)
    else
      pred <- predict(mymodel, newdata = outbootdata)
    if (!is.numeric(pred)) stop("predict does not return numerical values")
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")

    bootindx[-tindx, i] <- (pred - y[-tindx])^2
  }

  fun <- function(x)
        ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))

  err <- sqrt(mean(apply(bootindx, 1, fun), na.rm = TRUE))
  RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestreg"
  RET
}

bootest.Surv <- function(y, formula, data=NULL, model, predict, 
                           nboot=25, bc632plus=FALSE, ...) {
  
  # bootstrap estimator of Brier's score

  if (bc632plus) stop("cannot compute 632+ estimator of Brier's score")

  N <- dim(y)[1]
  nindx <- 1:N
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  USEPM <- FALSE

  if(is.null(nboot)) stop("number of bootstrap replications is missing")
  if (nboot <=2) stop("to small number of bootstrap replications")
  if (is.null(data)) data <- as.data.frame(rep(1, N))
  if (!is.data.frame(data)) stop("data is not a data.frame")

  for(i in 1:nboot) {
    tindx <- sample(nindx, N, replace = TRUE)
    mymodel <- model(formula, data=data[tindx,], ...)
    outbootdata <- subset(data, !(nindx %in% tindx))
    # check if mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM)
      pred <- predict(newdata=outbootdata)
    else
      pred <- predict(mymodel, newdata = outbootdata)

    if (is.list(pred)) {
      if (!inherits(pred[[1]], "survfit") && !inherits(pred, "survfit"))
        stop("predict does not return a list of survfit objects")
    } else {
      stop("predict does not return a list of survfit objects")
    }

    bootindx[-tindx, i] <- sbrier(y[-tindx], pred)
  }

  fun <- function(x)
        ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))

  err <- mean(apply(bootindx, 1, fun), na.rm = TRUE)
  RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestsurv"
  RET
}

