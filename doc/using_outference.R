## ------------------------------------------------------------------------
set.seed(2667)
n <- 100; p <- 6
# generate the design matrix
X <- matrix(rnorm(n*(p-1)), n, (p-1))
X <- scale(X, FALSE, TRUE)*sqrt(n/(n-1)) # col of X and has norm sqrt(n)
my.data <- as.data.frame(X)
X <- cbind(rep(1, n), X)
# generate the response with noise i.i.d. from N(0, 1)
beta <- rnorm(p)
y <- X %*% beta + rnorm(n)
# assume the first five observations are outliers, shifed upwards by 5 units
y[1:5] <- y[1:5] + 5
# format the dataset
my.data <- cbind(y = y, my.data)
head(my.data)

## ------------------------------------------------------------------------
library(outference)
fit <- outference(y ~ ., data = my.data, method = "cook", cutoff = 4, sigma = "estimate")
fit

## ---- fig.width = 6, fig.height = 4--------------------------------------
plot(fit)

## ------------------------------------------------------------------------
coef(fit)

## ------------------------------------------------------------------------
summary(fit)

## ------------------------------------------------------------------------
confint(fit, level = 0.95)

## ------------------------------------------------------------------------
# new data points for prediction
new.data <- t(rnorm(p-1))
colnames(new.data) <- colnames(my.data)[-1]
new.data <- as.data.frame(new.data)
predict(fit, newdata = new.data, interval = "confidence", level = 0.95)
predict(fit, newdata = new.data, interval = "prediction", level = 0.95)

## ------------------------------------------------------------------------
coeftest(fit, index = 2) # test the significance of "V1"
grptest(fit, group = 2:3) # test if "V1" and "V2" are simultaneously zero

## ------------------------------------------------------------------------
data("stackloss")
head(stackloss)

## ------------------------------------------------------------------------
require(outference)
stack.fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)

## ---- fig.width = 6, fig.height = 4--------------------------------------
plot(stack.fit)
# show which observations are considered as outliers
stack.fit$outlier.det

## ------------------------------------------------------------------------
summary(stack.fit)

## ------------------------------------------------------------------------
# full model
summary(stack.fit$fit.full)
# model with outliers removed, but inference is not corrected
summary(stack.fit$fit.rm)

## ------------------------------------------------------------------------
confint(stack.fit, level = 0.95)

## ------------------------------------------------------------------------
confint(stack.fit$fit.rm, level = 0.95)

## ------------------------------------------------------------------------
# say we want to predict with the first observation 
newdata = stackloss[1, 1:3]
# CI for regression surface (a.k.a., short prediction intervals)
predict(stack.fit, newdata = newdata, interval = "confidence", level = 0.95)
# prediction intervals (a.k.a., wide prediction intervals) 
predict(stack.fit, newdata = newdata, interval = "prediction", level = 0.95)

# predict with the whole dataset is also possible by running the following command, 
# but it may take a while

# predict(stack.fit, interval = "confidence", level = 0.95)

## ---- fig.width = 6, fig.height = 4--------------------------------------
stack.lasso = outference(stack.loss ~ ., data = stackloss, method = "lasso")
plot(stack.lasso)
stack.lasso$outlier.det

## ------------------------------------------------------------------------
summary(stack.lasso)

## ---- fig.width = 6, fig.height = 4--------------------------------------
stack.dffits = outference(stack.loss ~ ., data = stackloss, method = "dffits")
plot(stack.dffits)
stack.dffits$outlier.det

## ------------------------------------------------------------------------
summary(stack.dffits)

