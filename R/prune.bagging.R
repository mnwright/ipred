# $Id: prune.bagging.R,v 1.1 2002/03/02 14:21:15 hothorn Exp $

prune.bagging <- function(tree, cp=0.01,...)
{
  for(i in 1:length(tree$mt))
    tree$mt[i][[1]] <- prune(tree$mt[i][[1]], cp=cp, ...)
  tree
}
