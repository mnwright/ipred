# $Id: zzz.R,v 1.8 2002/09/12 08:59:13 hothorn Exp $

.First.lib <- function(lib, pkg) {
    if(!require(rpart))
        warning("Could not load package rpart")
    if(!require(MASS))
        warning("Could not load package MASS")
    if(!require(mlbench))
        warning("Could not load package mlbench")
    if(!require(survival))
        warning("Could not load package mlbench")
    if(!require(class))
        warning("Could not load package class")
    if(!require(nnet))
        warning("Could not load package nnet")
}
