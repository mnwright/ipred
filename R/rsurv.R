# $Id: rsurv.R,v 1.3 2003/02/13 10:49:41 hothorn Exp $

rsurv <- function(N, model=c("A", "B", "C", "D", "tree"), gamma=NULL, fact=1,
                  pnon=10)
{
    model <- match.arg(model)
    X <- matrix(runif(N*5), ncol=5)
    colnames(X) <- paste("X", 1:ncol(X), sep="")
    switch(model,
        "A" =  { 
            time <- rexp(N)
        }, 
        "B" = {
            hazard <- as.numeric(X[,1] <= 0.5 & X[,2] > 0.5)
            time <- rexp(N)
            time[hazard == 1] <- rexp(sum(hazard==1), exp(3))
        },
        "C" = {
            hazard <- 3*X[,1] + X[,2]
            time <- sapply(exp(hazard), rexp, n=1)
        },
        "D" = {
            hazard <- 3*X[,1] - 3*X[,2] + 4*X[,3] - 2*X[,4]
            time <- sapply(exp(hazard), rexp, n=1)
        },
        "tree" = {
            hazard <- rep(0, nrow(X))
            hazard[(X[,1] <= 0.5 & X[,2] <= 0.5)] <- 0
            hazard[(X[,1] <= 0.5 & X[,2] > 0.5 & X[,4] <= 0.5)] <- 1
            hazard[(X[,1] <= 0.5 & X[,2] > 0.5 & X[,4] > 0.5)] <- 0
            hazard[(X[,1] > 0.5 & X[,3] <= 0.5 & X[,5] <= 0.3)] <- 1
            hazard[(X[,1] > 0.5 & X[,3] <= 0.5 & X[,5] > 0.3)] <- 2
            hazard[(X[,1] > 0.5 & X[,3] > 0.5 & X[,4] <= 0.7)] <- 2
            hazard[(X[,1] > 0.5 & X[,3] > 0.5 & X[,4] > 0.7)] <- 3
            hazard <- hazard * fact
            time <- sapply(exp(hazard), rexp, n=1)
            if (pnon > 0)
              X <- cbind(X, matrix(runif(N*pnon), ncol=pnon))
            colnames(X) <- paste("X", 1:ncol(X), sep="")
        })
    if (!is.null(gamma))  
        censtime <- runif(N, min=0, max=gamma)
    else
        censtime <- Inf
    cens <- as.numeric(time <= censtime)
    time <- pmin(time, censtime)
    simdat <- as.data.frame(cbind(time, cens, X))
    return(simdat)
}
    