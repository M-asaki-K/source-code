#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(11:14)]

#-----------remove columns of 0 distribution from x----
x.sds <- apply(x, 2, sd)
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]

#-----------standarization of y------------------------
preprocessed.y <- (y - mean(y)) / sd(y)
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x------------------------
apply(x, 2, mean)
apply(x, 2, sd)
preprocessed.x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
preprocessed.x <- data.frame(preprocessed.x)

#-----------compare the number of columns and rows--------------
ncol(preprocessed.x)
nrow(preprocessed.x)

#-----------pick up columns if needed---------------------------
multi.regression.x <- preprocessed.x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(preprocessed.y, multi.regression.x)

#-----------verify each factor's relations--------------------------
multi.regression.compounds <- as.data.frame(multi.regression.compounds)

#-----------MLR regression--------------------------------------------
compounds.lm <- lm(preprocessed.y~., data=multi.regression.compounds)
compounds.lm

summary(compounds.lm)

lm.predicted.y <- predict(compounds.lm)
lm.predicted.y

cor(preprocessed.y, lm.predicted.y)
lm.r2 <- cor(preprocessed.y, lm.predicted.y)**2
lm.r2

plot(preprocessed.y, lm.predicted.y,
     xlab="Observed value",
     ylab="Predicted value")
abline(a=0, b=1)

res <- lm.predicted.y - preprocessed.y
y2 <- res
preprocessed.y2 <- (y2 - mean(y2)) / sd(y2)

x2 <- complete.compounds[,-c(1,11:14)]
preprocessed.x2 <- apply(x2, 2, function(x2) {(x2 - mean(x2)) / sd(x2)})
preprocessed.x2 <- data.frame(preprocessed.x2)

#-----------pick up columns if needed---------------------------
multi.regression.x2 <- preprocessed.x2[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds2 <- cbind(preprocessed.y2, multi.regression.x2)

#-----------verify each factor's relations--------------------------
multi.regression.compounds2 <- as.data.frame(multi.regression.compounds2)

#-----------MLR regression--------------------------------------------
compounds.lm2 <- lm(preprocessed.y2~., data=multi.regression.compounds2)
compounds.lm2

summary(compounds.lm2)

lm.predicted.y2 <- predict(compounds.lm2)
lm.predicted.y2

cor(preprocessed.y2, lm.predicted.y2)
lm.r2 <- cor(preprocessed.y2, lm.predicted.y2)**2
lm.r2

plot(preprocessed.y2, lm.predicted.y2,
     xlab="Observed value",
     ylab="Predicted value")
abline(a=0, b=1)

#-----------definition of scaled.compounds (used for PLS)-----------
scaled.compounds2 <- cbind(preprocessed.y2, preprocessed.x2)
scaled.compounds2 <- as.data.frame(scaled.compounds2)

#-------ここまでデータ前処理--------
#PLS
library(pls)

compounds.plsr2 <- plsr(preprocessed.y2~., data=scaled.compounds2, validation="CV")
summary(compounds.plsr2)
plot(compounds.plsr2, "validation")
compounds.plsr$validation$PRESS[1,]

ncomp.onesigma2 <- selectNcomp(compounds.plsr2, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma2

predict(compounds.plsr2)[, , ncomp.onesigma2]
plsr.predicted.y2 <- predict(compounds.plsr2)[, , ncomp.onesigma2]
plsr.r22 <- cor(preprocessed.y2, plsr.predicted.y2)**2
plsr.r22

plot(preprocessed.y2, plsr.predicted.y2,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLSR")
abline(a=0, b=1)

compounds.plsr2$coefficients[, , ncomp.onesigma2]


library(plsVarSel)

#----------ここからVIP-------------------------
comp <- which.min(compounds.plsr2$validation$PRESS)
vip <- VIP(compounds.plsr2,comp)
vip
plot(vip,compounds.plsr2$coefficients[, , ncomp.onesigma2])
vip.selected <- bve_pls(preprocessed.y2, preprocessed.x2, ncomp = ncomp.onesigma2, VIP.threshold = 1)
#usually ncomp = 10 as default, but in that case ncomp is optimized by PRESS (predicted residual sum of squares), 
#which tends to overcount the optimum number of components.
#thus here, maximum of ncomp is determined exceptionally by randomization method.
vip.selected

x3 <- complete.compounds[,c(2,10,15)]
View(x3)

#-----------standarization of y------------------------
y3 <- complete.compounds[,c(1)]
y3 <- (y3 - mean(y3)) / sd(y3)
preprocessed.y3 <- (y3 - min(y3)) / (max(y3) - min(y3))

#-----------standarization of x------------------------
x3 <- apply(x3, 2, function(x3) {(x3 - mean(x3)) / sd(x3)})
preprocessed.x3 <- apply(x3, 2, function(x3) {(x3 - min(x3)) / (max(x3) - min(x3))})
View(preprocessed.x3)

#-----------x converted into data frame type for machine learning-----------
preprocessed.x0 <- data.frame(preprocessed.x3[,])
View(preprocessed.x0)

preprocessed.x1 <- preprocessed.x0[c(1:700),]
preprocessed.x2 <- preprocessed.x0[c(701:800),]
preprocessed.y1 <- preprocessed.y3[1:700]
preprocessed.y2 <- preprocessed.y3[701:800]
View(preprocessed.x2)

library(GPfit)
x = preprocessed.x1;
View(preprocessed.x1)
y = preprocessed.y1;
GPmodel = GP_fit(x,y);
print(GPmodel)

GPprediction = predict.GP(GPmodel,preprocessed.x2);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
yhat

plot(yhat,preprocessed.y2)

sigma <- 3*(mse**0.5)

f <- cbind(yhat, preprocessed.y2, sigma)
f <- as.data.frame(f)

library(reshape2)
library(ggplot2)
library(ggsci)

g <- ggplot(f, aes(x = preprocessed.y2, y = yhat))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin = yhat - sigma, ymax = yhat + sigma, width = 0.03))
g <- g + scale_fill_nejm()
g <- g + geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1)

plot(g)

#Partial dependence plot

m <- c(-0.5,-0.2,0.2,0.4,0.6,0.8,1.2,1.5) # valiable
n <- c(0,0,0,0,0,0,0,0) # constant
M.data <-cbind.data.frame(n,n,m) #you can choose the valiable manually

GPprediction = predict.GP(GPmodel,xnew=M.data);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
plot(m,yhat)

sigma <- 3*(mse**0.5) # 3 sigma
f <- cbind(yhat, m, sigma)
f <- as.data.frame(f)

g <- ggplot(f, aes(x = m, y = yhat))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin = yhat - sigma, ymax = yhat + sigma, width = 0.03))
g <- g + scale_fill_nejm()
g <- g + xlim(-1,2) + ylim(0.4,1.4)

plot(g)

#----------------------------Bayesian optimization_grid search--------------------------------------
parameter_list <- list(a1 = seq(0,1,length.out = 10),
                       a2 = seq(0,1,length.out = 10),
                       a3 = seq(0,1,length.out = 10)
                      
)

grid <- unname(expand.grid(parameter_list, KEEP.OUT.ATTRS = FALSE)) #form grid
View(grid)
pred <- predict.GP(GPmodel,grid)
acq_func <- function(mu, var, k) {
  acqs <- mu - k * sqrt(var)
  which.min(acqs)
}

ind <-grid[acq_func(pred$Y_hat,pred$MSE,3),]
ind #print optimum parameter
min(pred$Y_hat-3*pred$MSE**0.5) #print maximum

#---------------------Bayesian optimization_find optim-------------
fn <- function(a){
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]

  model <- GPmodel
  t.Pred <- predict(model, data.frame(a1, a2, a3))
  (t.Pred$Y_hat - 3*(t.Pred$MSE)**0.5)}

optim(c(1,1,0), fn, method = "L-BFGS-B", lower = c(0,0,0), upper = c(1,1,1),control=list(fnscale=1))

#---------------------------optimization------------------------
library(BB)
spg(c(1,1,0), fn, method = 3, lower = c(0,0,0), upper = c(1,1,1),control=list(maximize=FALSE))

