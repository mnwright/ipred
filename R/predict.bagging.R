# $Id: predict.bagging.R,v 1.7 2002/05/16 15:22:11 hothorn Exp $

uwhich.max <- function(x) {
  wm <- which.max(x)
  if (length(wm) > 1)
    wm <- wm[sample(length(wm), 1)]
  wm
}

predict.bagging <- function(object, newdata=NULL, ...) 
{
    if (!inherits(object, "bagging")) {
        if(inherits(object,"rpart")) {
            if (object$method == "class")
                RET <- predict.rpart(object, newdata, type="class") 
            else
                RET <- predict.rpart(object, newdata)
        }
    } else {
        LDA <- FALSE
        if (is.null(newdata)) stop("newdata missing")
        if (!is.null(object$mt)) mt <- object$mt
        if (!is.null(object$ldasc)) { LDA <- TRUE; ldasc <- object$ldasc; }
        if (mt[[1]]$method == "class") {
            classlevels <- attr(mt[[1]], "ylevels")
            votenew <- matrix(0, nrow=nrow(newdata),
                                 ncol=length(classlevels))
            for (i in 1:length(mt)) {
                if (LDA) { 
		    ldap <- rpart.matrix(newdata)
                    class(ldap) <- NULL
                    test <- cbind(newdata, predict(ldasc[[i]],
                                           newdata=as.data.frame(ldap))$x)
                } else {
                    test <- newdata
                }
                pr <- predict.rpart(mt[[i]], test, type="class")
                votenew[cbind(1:nrow(votenew), as.integer(pr))] <-
                    votenew[cbind(1:nrow(votenew), as.integer(pr))] + 1
            }
            RET <- apply(votenew, 1, uwhich.max)
            RET <- as.factor(RET)
            levels(RET) <- classlevels
        } else {
            if (!is.null(object$mt)) mt <- object$mt
            if (!is.null(object$ldasc)) 
              stop("cannot predict with lda for regression trees!")
            a <- predict.rpart(mt[[1]], newdata)
            for (i in 2:length(mt))
                a <- a + predict.rpart(mt[[i]], newdata)
            RET <- a/length(mt)
        }
    }  
    RET
}
