# $Id: predict.inclass.R,v 1.15 2002/04/05 14:36:09 peters Exp $

# Additional option type ="class", if intermediate is nominal

predict.inclass <- function(object, cFUN, intbag = NULL, newdata, ...)
{
  q <- length(object)		# number of intermediates
  classes <- c()
  for(i in 1:q) {classes <- c(classes, class(object[[i]]))}
  result <- c()
  res <- c()

  if(!is.null(intbag) && !all(classes == "bagging")) {
    stop("intbag is only specified for intermediates of class bagging")
  }			

  if(!is.null(intbag) && intbag == FALSE && all(classes == "bagging")) {			
    # first classification, than bagging
    diagbag <- list()
    # Calculation of a matrix with estmations for each bootstrap sample and 
    # each intermediate variable 
    for(i in 1:q) {		# over intermediates
      if(class(object[[i]])=="bagging") object[[i]] <- object[[i]]$mt  
      K <- length(object[[i]])	# number of bootstrap samples
      dummy <- list()  
      for(j in 1:K) {		# over bagging
        if(inherits(object[[i]][[j]], "rpart")) {
          if(object[[i]][[j]]$method == "class")
            RET <- predict.rpart(object[[i]][[j]], newdata, type="class")
          else
            RET <- predict.rpart(object[[i]][[j]], newdata)
        }
        dummy <- c(dummy, list(RET))        
      }
      diagbag <- c(diagbag,list(dummy))
    }
  diagbag <- as.data.frame(diagbag)
  # Calculation of the diagnosis for each bootstrap sample   
    for(i in 1:K) {
      dummy <- diagbag[,i+((0:(q-1))*K)]
      names(dummy) <- names(object)
      res <- cbind(res, as.character(cFUN(data = as.data.frame(dummy))))
    }
    # Voting
    vfun <- function(votes) {
      votes <- as.factor(votes)
      levels(votes)[which.max(table(votes))]
    }
    result <- as.factor(apply(as.data.frame(res), 1, vfun))
  }

  if(intbag == TRUE || is.null(intbag)) {			
    # first bagging, than classification
    intermediate <- list()
    for(i in 1:q) {
      if(inherits(object[[i]], "lm")) {
        RET <- mypredict.lm(object[[i]], newdata = newdata)
      }
      if(inherits(object[[i]], "rpart")) {
        if(object[[i]]$method == "class")
          RET <- predict.rpart(object[[i]], newdata, type="class")
        else
          RET <- predict.rpart(object[[i]], newdata)
      }
      if(inherits(object[[i]], "bagging")) {
        RET <- predict(object[[i]], newdata = newdata, ...)
      }
      if(inherits(object[[i]], "lda")) {
        RET <- predict(object[[i]], newdata, ...)$class
      }
      intermediate <- c(intermediate, list(RET))
    }

    names(intermediate) <- names(object)
    intermediate <- as.data.frame(intermediate)
    result <- cFUN(intermediate)
  }
  return(result)
}
