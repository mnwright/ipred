# $Id: zzz.R,v 1.12 2003/07/22 14:56:31 peters Exp $

.onLoad <- function(lib, pkg) {
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
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
    library.dynam("ipred", pkg, lib)
}
