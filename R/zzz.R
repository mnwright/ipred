# $Id: zzz.R,v 1.7 2002/03/26 16:29:15 hothorn Exp $

.First.lib <- function(lib, pkg) {
    if(!require(rpart))
        warning("Could not load package rpart")
    if(!require(MASS))
        warning("Could not load package MASS")
    if(!require(mlbench))
        warning("Could not load package mlbench")
}
