# $Id: predict.inclass.R,v 1.18 2003/03/12 17:06:25 hothorn Exp $

# Additional option type ="class", if intermediate is nominal

predict.inclass <- function(object, cFUN, intbag = NULL, newdata, ...)
{
  q <- length(object)		# number of intermediates
  classes <- c()
  for(i in 1:q) {classes <- c(classes, class(object[[i]]))}
  result <- c()
  res <- c()

  if(!is.null(intbag) && !all(classes == "regbagg" | classes == "classbagg")) {
    stop("intbag is only specified for intermediates of class bagging")
  }			

  if(!is.null(intbag) && intbag == FALSE && all(classes == "regbagg" | classes == "classbagg")) {			
    # first classification, than bagging
    diagbag <- list()
    # Calculation of a matrix with estmations for each bootstrap sample and 
    # each intermediate variable 
    for(i in 1:q) {		# over intermediates
      if(classes == "regbagg" | classes == "classbagg") 
        object[[i]] <- object[[i]]$mtrees  
      K <- length(object[[i]])	# number of bootstrap samples
      dummy <- list()  
      for(j in 1:K) {		# over bagging
        if(inherits(object[[i]][[j]], "rpart")) {
          if(object[[i]][[j]]$method == "class")
            RET <- predict(object[[i]][[j]], newdata, type="class")
          else
            RET <- predict(object[[i]][[j]], newdata)
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
          RET <- predict(object[[i]], newdata, type="class")
        else
          RET <- predict(object[[i]], newdata)
      }
      if(inherits(object[[i]],  "regbagg") | inherits(object[[i]], "classbagg")) {
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
