library(ipred)
data(BreastCancer)
set.seed(29081975)

mod <- bagging(Class ~ Cl.thickness + Cell.size
                + Cell.shape + Marg.adhesion
                + Epith.c.size + Bare.nuclei
                + Bl.cromatin + Normal.nucleoli
                + Mitoses, data=BreastCancer, coob=TRUE)
print(mod)

print(predict(mod, newdata=BreastCancer))

# bagging failed if only one predictor was specified
# by Christoph M. Friedrich <chris@uni-wh.de>, April 29th, 2002

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))
learn <- cbind(y, X)
mt <- bagging(y ~ V1, data=learn)
mt <- bagging(y ~ V1, data=learn, method="double", coob=FALSE)
X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- apply(X, 1, mean) + rnorm(nrow(X))
learn <- cbind(y, X)
mt <- bagging(y ~ V1, data=learn)

