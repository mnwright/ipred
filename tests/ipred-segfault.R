library('ipred')
set.seed(29081975)

# Classification

learn <- as.data.frame(mlbench.twonorm(200))
test <- as.data.frame(mlbench.twonorm(100))

# bagging

mod <- bagging(classes ~ ., data=learn, coob=TRUE, nbagg=10)
mod
predict(mod)[1:10]

# Double-Bagging

comb.lda <- list(list(model=lda, predict=function(obj, newdata)
                      predict.lda(obj, newdata)$x))

mod <- bagging(classes ~ ., data=learn, comb=comb.lda, nbagg=10)
mod
predict(mod, newdata=test[1:10,])
predict(mod, newdata=test[1:10,], agg="aver")
predict(mod, newdata=test[1:10,], agg="wei")
predict(mod, newdata=test[1:10,], type="prob")
predict(mod, newdata=test[1:10,], type="prob", agg="aver")
predict(mod, newdata=test[1:10,], type="prob", agg="wei")

mypredict.lda <- function(object, newdata)
       predict(object, newdata = newdata)$class

errorest(classes ~ ., data=learn, model=lda, predict=mypredict.lda)
errorest(classes ~ ., data=learn, model=lda, predict=mypredict.lda,
est.para=list(k=5, random=FALSE))
errorest(classes ~ ., data=learn, model=bagging, est.para=list(k=2),
nbagg=10)
errorest(classes ~ ., data=learn, model=bagging, est.para=list(k=2),
nbagg=10, comb=comb.lda)
errorest(classes ~ ., data=learn, model=lda,
predict=mypredict.lda, estimator="boot")
errorest(classes ~ ., data=learn, model=lda,
predict=mypredict.lda, estimator="632plus")

# Regression

learn <- as.data.frame(mlbench.friedman1(100))
test <- as.data.frame(mlbench.friedman1(100))

# bagging

mod <- bagging(y ~ ., data=learn, coob=TRUE, nbagg=10)
mod
predict(mod)[1:10]

predict(mod, newdata=test[1:10,])
predict(mod, newdata=test[1:10,], agg="aver") 
predict(mod, newdata=test[1:10,], agg="wei")  
errorest(y ~ ., data=learn, model=lm)
errorest(y ~ ., data=learn, model=lm, est.para=list(k=5, random=FALSE))
errorest(y ~ ., data=learn, model=lm, estimator="boot")

# survival

learn <- rsurv(100, model="C")
test <- rsurv(100, model="C")

mod <- bagging(Surv(time, cens) ~ ., data=learn, nbagg=10)
mod
predict(mod, newdata=test[1:10,])

errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
         est.para=list(k=2, random=FALSE), nbagg=5)
errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
         estimator="boot", nbagg=5, est.para=list(nboot=5))

# bundling for regression

learn <- as.data.frame(mlbench.friedman1(100))
test <- as.data.frame(mlbench.friedman1(100))

comb <- list(list(model=lm, predict=predict.lm))

modc <- bagging(y ~ ., data=learn, nbagg=10, comb=comb)
modc
predict(modc, newdata=learn)[1:10]

# bundling for survival

data(GBSG2)
rcomb <- list(list(model=coxph, predict=predict.coxph))

mods <- bagging(Surv(time,cens) ~ ., data=GBSG2, nbagg=10, 
                comb=rcomb,  control=rpart.control(xval=0))
predict(mods, newdata=GBSG2[1:3,])
