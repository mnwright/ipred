# $Id: bagging.R,v 1.6 2002/03/28 07:42:25 peters Exp $

bagging <- function(y, ...) UseMethod("bagging")

bagging.default <- function(y, X=NULL, nbagg=25, method=c("standard","double"),
                            coob=TRUE, control=rpart.control(minsize=2, cp=0),
                            ...) {
  method <- match.arg(method)
  class <- is.factor(y)
  if (method=="double" & !class)
      stop("Cannot compute Double Bagging for regression problems")
  if (!is.data.frame(X)) stop("X is not a data.frame")
  mt <- list()
  ldasc <- c()
  oob <- list()
  for (i in 1:nbagg) {
    indx <- sample(1:length(y), length(y), replace=TRUE)
    if (method == "standard")
      mt <- c(mt, list(rpart(y[indx] ~., data=X[indx,], control = control,...)))
    if (method == "double") {
      lmod <- lda(y[-indx] ~., data=X[-indx,]) # OOB !!!
      ldasc <- c(ldasc, list(lmod$scaling))
      mt <- c(mt, list(rpart(y[indx] ~.,
                  data=cbind(as.matrix(X[indx,])%*%lmod$scaling, X[indx,]),
                  control=control, ...)))
    }
    if (coob & method != "double") {
      help <- predict.bagging(mt[[i]], X[-indx,])
      if (any(is.na(help))) warning("NA in predict")
      dummy <- rep(NA, length(y))
      if (class) {
        dummy <- as.factor(dummy)
        levels(dummy) <- levels(help)
      }
      dummy[-indx] <- help
      oob <- c(oob, list(dummy))
    }
  }
  if (coob & method == "double") {
        warning("Cannot compute out-of-bag estimate for Double-Bagging")
        err <- NA
        oob <- NA
  }
  if (coob & method =="standard") {
    oob <- as.data.frame(oob)
    names(oob) <- paste("t", 1:length(mt), sep="")
    vfun <- function(votes) {
      if (all(is.na(votes))) return(NA)
      votes <- votes[!is.na(votes)]
      votes <- as.factor(votes)
      levels(votes)[which.max(table(votes))]
    }
    if (class) {
        pred <- as.factor(unlist(apply(oob, 1, vfun)))
        err <- mean(y != pred, na.rm=TRUE)
    } else {
        pred <- apply(as.data.frame(oob), 1, mean, na.rm=TRUE)
        err <- mean( (y - pred)^2, na.rm=TRUE)
    }
  } else { 
    err <- NA
    oob <- NA
  }
  if (method =="double")
    method <- "Double-Bagging"
  else
    method <- "Bagging"
  if (class)
    method <- paste(method, "classification trees")
  else 
    method <- paste(method, "regression trees")
  RET <- list(mt=mt, oob=oob, err=err, nbagg=nbagg, method=method)
  if (method == "double") RET <- list(mt=mt, oob=oob, err=err, ldasc = ldasc,
                                      nbagg=nbagg, method=method)
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
    DATA <- list(y = mf[[response]], X = as.data.frame(X))
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
