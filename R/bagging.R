# $Id: bagging.R,v 1.14 2002/06/17 08:53:00 hothorn Exp $

bagging <- function(y, ...) UseMethod("bagging")

bagging.default <- function(y, X=NULL, nbagg=25, method=c("standard","double"),
                            coob=FALSE, control=rpart.control(minsplit=2, cp=0),
                            tol=1.0e-4, ...) {
  method <- match.arg(method)
  class <- is.factor(y)
  if (class & coob) {
    classlevels <- levels(y)
    votenew <- matrix(0, nrow=length(y), ncol=length(classlevels))
  }
  oobsum <- 0
  if (method=="double" & !class)
      stop("Cannot compute Double Bagging for regression problems")
  if (!is.data.frame(X)) stop("X is not a data.frame")
  mt <- list()
  ldasc <- c()
  for (i in 1:nbagg) {
    indx <- sample(1:length(y), length(y), replace=TRUE)
    if (method == "standard")
      mt <- c(mt, list(rpart(y[indx] ~., data=as.data.frame(X[indx,]),
                             control = control,...)))
    if (method == "double") {
      # there are some problems with predict.lda for only one predictor.
      if (ncol(X) < 2) 
        stop("cannot compute double-bagging with less than two predictors")
      lday <- y[-indx]
      ldap <- as.data.frame(X[-indx,])
      lmod <- lda(lday ~ ., data=ldap, tol=tol) # OOB !!!
      ldasc <- c(ldasc, list(lmod))
      ldap <- as.data.frame(X[indx,])
      rpartp <- cbind(ldap, predict(lmod, newdata=ldap)$x)
      mt <- c(mt, list(rpart(y[indx] ~ ., data=as.data.frame(rpartp), 
                             control=control, ...)))
    }
    if (coob & method != "double") {
      if (class) {
        pr <- predict.bagging(mt[[i]], newdata=X, type="class")
        votenew[cbind((1:length(y))[-indx], as.integer(pr[-indx]))] <-
          votenew[cbind((1:length(y))[-indx], as.integer(pr[-indx]))] + 1
      } else {
        pr <- predict.bagging(mt[[i]], newdata=X)
        oobsum <- oobsum + pr 
      }
      if (any(is.na(pr))) warning("NA in predict")
    }
  }
  if (coob & method == "double") {
        warning("Cannot compute out-of-bag estimate for Double-Bagging")
        err <- NA
        pred <- NA
  }
  if (coob & method =="standard") {
    if (class) {
        pred <- apply(votenew, 1, uwhich.max) 
        pred <- as.factor(pred)
	levels(pred) <- classlevels
        err <- mean(y != pred, na.rm=TRUE)
    } else {
        pred <- oobsum/nbagg
        err <- mean( (y - pred)^2, na.rm=TRUE)
    }
  } else { 
    err <- NA
    pred <- NA
  }
  if (method =="double")
    mmethod <- "Double-Bagging"
  else
    mmethod <- "Bagging"
  if (class)
    mmethod <- paste(mmethod, "classification trees")
  else 
    mmethod <- paste(mmethod, "regression trees")
  RET <- list(mt=mt, oob=pred, err=err, nbagg=nbagg, method=mmethod)
  if (method == "double") RET <- list(mt=mt, oob=pred, err=err, ldasc = ldasc,
                                      nbagg=nbagg, method=mmethod)
  class(RET) <- "bagging"

  return(RET)
}


bagging.formula <-
function(formula, data, subset, na.action=na.rpart, ...)
{
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) < 1)
       || (length(attr(terms(formula[-3]), "term.labels")) != 1))
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    if (length(names(mf)) > 2) {
        DNAME <- paste(names(mf)[1], "by", paste(names(mf)[-1],
                 collapse = " + "))
    } else {
        DNAME <- paste(names(mf), collapse = " by ")
    }
    response <- attr(attr(mf, "terms"), "response")
    X <- rpart.matrix(mf)
    class(X) <- NULL
    DATA <- list(y = mf[[response]], X = as.data.frame(X)) # mf[,-1]))
    names(DATA) <- c("y", "X")
    y <- do.call("bagging", c(DATA, list(...)))
    y$data.name <- DNAME
    return(y)
}

print.bagging <- function(x, digits=4, ...)
{
    cat("\n")
    method <- x$method 
    method <- paste(method, "\nwith", x$nbagg, "bootstrap replications")
    cat(method, "\n")
    cat("\n")
    cat("Data: ", x$data.name, "\n")
    cat("\n")
    if (x$mt[1][[1]]$method == "class")
        error <- "Out-of-bag misclassification error:"
    else
        error <- "Out-of-bag mean squared error:"
    if (!is.na(x$err))
        cat(error, round(x$err, digits), "\n")
}

summary.bagging <- function(object, ...)
{
     class(object) <- "summary.bagging"
     object
}

print.summary.bagging <- function(x, digits = max(3, getOption("digits")-3),
                                 ...)
{
     print.bagging(x, digits=digits, ...)
     cat("\n")
     cat("Trees: \n")
     print(x$mt)
     invisible(x$mt)
}
