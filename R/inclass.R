# $Id: inclass.R,v 1.24 2002/09/24 15:22:20 peters Exp $

####################################################
# options for predict.inclass are only: rpart, lm, bagging
# Assumption: All intermediate variables have the same level

inclass <- function(object, ...) UseMethod("inclass")

inclass.default <- function(object, ...)
{
  stop(paste("Do not know how to handle objects of class", class(object)))
}



flist <- function(formula, ...){
 if(missing(formula)) stop("formula missing")

 m <- match.call(expand.dots = FALSE)

 if(is.list(formula) && is.null(m$...)) {
   result <- formula
 } 

 if(is.list(formula) && !is.null(m$...)) {
   result <- c(eval(m$formula), m$...)
 }

 if(!is.list(formula)) {
   if((length(formula) != 3)
     || (length(attr(terms(formula[-2]), "term.labels")) < 1))
   stop("formula missing or incorrect")

   q <- length(attr(terms(formula[-3]), "term.labels"))
   formula.list <- list()
   for(i in 1:q) {
    formelnew <- formula
    formelnew[[2]] <- as.name(attr(terms(formula[-3]), "term.labels")[i])
    formula.list <- c(formula.list, formelnew)
  }
  result <- c(formula.list, ...)
 }
 
 class(result) <- "flist"
 return(result)
 }

inclass.formula <- function(formula, pFUN, data, subset, na.action, ...) 
{
  if(missing(formula)
    || (length(formula) != 3)
    || (length(attr(terms(formula[-2]), "term.labels")) < 1))
    stop("formula missing or incorrect")

  q <- length(attr(terms(formula[-3]), "term.labels"))
  formula.list <- list()
  for(i in 1:q) {
    formelnew <- formula
    formelnew[[2]] <- as.name(attr(terms(formula[-3]), "term.labels")[i])
    formula.list <- c(formula.list, formelnew)
  }

  rm(formelnew, formula)
  formula.list <- flist(formula.list)
  result <- inclass.flist(object = formula.list, pFUN=pFUN, data = data, subset = subset, na.action = na.action, ...)
  return(result)
}


inclass.flist <- function(object, pFUN, data, subset, na.action, ...) 
{
  formula.list <- object
  q <- length(formula.list)
  responses <- c()
  for(i in 1:q) {
    responses <- c(responses, paste(formula.list[[i]][[2]]))
  }

  result <- list()
  namen <- c()

  for(i in 1:q) {
    formula <- formula.list[[i]]
    data.without <- data[,!(names(data) %in% responses[-i])]
    if(missing(formula)
      || (length(formula) != 3)
      || (length(attr(terms(formula[-2]), "term.labels")) < 1)
      || (length(attr(terms(formula[-3]), "term.labels")) != 1))
      stop("formula missing or incorrect")
    m <- match.call(expand.dots= FALSE)
    if(missing(subset)) {
      if(missing(na.action))
        res <- pFUN(formula = formula, data = data.without, ...)
      else 
        res <- pFUN(formula = formula, data = data.without, na.action = na.action, ...)
    } else {
      if(missing(na.action))
        res <- pFUN(formula = formula, data = data.without, subset = subset, ...)
      else
        res <- pFUN(formula = formula, data = data.without, subset = subset, na.action = na.action, ...)
    }
    namen <- c(namen, as.character(formula[[2]]))
    result <- c(result, list(res))
  }

  names(result) <- namen
  class(result) <- "inclass"
  return(result)
}


print.inclass <- function(x, ...)
{
  q <- length(x)
  intermediates <- attr(x, "names")
  classes <- class(x[[1]])
 
  text.intermediates <- paste("Indirect classification, with", q, "intermediate variables:")
  predictive  <- paste("Predictive model per intermediate is", classes)
  predictive <- ifelse(classes == "bagging", paste(predictive, "with", x[[1]]$nbagg, "bootstrap replications"), predictive)
 
  cat("\n", text.intermediates, "\n", intermediates, "\n", "\n", predictive, "\n") 
}


summary.inclass <- function(object, ...)
{
  class(object) <- "summary.inclass"
  object
}


print.summary.inclass <- function(x, ...)
{
  q <- length(x)
  intermediates <- attr(x, "names")
  classes <- class(x[[1]])
 
  text.intermediates <- paste("Indirect classification, with", q, "intermediate variables:")
  predictive  <- paste("Predictive model per intermediate is", classes)
  predictive <- ifelse(classes == "bagging", paste(predictive, "with", x[[1]]$nbagg, "bootstrap replications"), predictive)

  cat("\n", text.intermediates, "\n", intermediates, "\n", "\n", predictive, "\n", "\n", "Models:", "\n")
  for(i in 1:length(x)) { 
    print(x[[i]])
  }
}
 