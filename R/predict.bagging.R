# $Id: predict.bagging.R,v 1.4 2002/03/28 07:42:25 peters Exp $

predict.bagging <- function(object, newdata=NULL, ...) 
{
    a <- c()
    if (!inherits(object, "bagging")) {
        if(inherits(object,"rpart")) {
            if (object$method == "class")
                RET <- predict.rpart(object, newdata, type="class") 
            else
                RET <- predict.rpart(object, newdata)
        }
    } else {
        LDA <- FALSE
        if (!is.null(object$mt)) mt <- object$mt
        if (!is.null(object$ldasc)) { LDA <- TRUE; ldasc <- object$ldasc; }
        if (mt[1][[1]]$method == "class") {
            a <- list()
            for (i in 1:length(mt)) {
                if (LDA) { 
                    test <- cbind(as.matrix(newdata)%*%ldasc[i][[1]], newdata)
                    names(test)[1] <- "LD1"
                } else {
                    test <- newdata
                }
                a <- c(a,list(predict.rpart(mt[i][[1]], test, type="class")))
            }
            names(a) <- paste("t", 1:length(mt), sep="")
            vfun <- function(votes) {
                votes <- as.factor(votes)
                levels(votes)[which.max(table(votes))]
            }
            RET <- as.factor(apply(as.data.frame(a), 1, vfun))
        } else {
            if (!is.null(object$mt)) mt <- object$mt
            if (!is.null(object$ldasc)) 
              stop("cannot predict with lda for regression trees!")
            for (i in 1:length(mt))
                a <- cbind(a,predict.rpart(mt[i][[1]], newdata))
            RET <- apply(a, 1, mean)
        }
    }  
    RET
}

