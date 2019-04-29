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
x.0 <- complete.compounds[,-c(1)]
x.2 <- apply(x.0,2,function(x.0) {x.0**3}) #add nonlinear columns
View(x.2)

x <- cbind.data.frame(x.0,x.2)


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

#--------------------divide into test and training data----------------------
train_size = 0.1

n = nrow(multi.regression.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[perm, ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-perm, ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)

#----------------------MLR regression training--------------------------------------------
compounds.lm <- lm(preprocessed.y~., data=multi.regression.compounds.train)
compounds.lm

summary(compounds.lm)

lm.predicted.y <- predict(compounds.lm)
lm.predicted.y

cor(preprocessed.y.train, lm.predicted.y)
lm.r2 <- cor(preprocessed.y.train, lm.predicted.y)**2
lm.r2

plot(preprocessed.y.train, lm.predicted.y,
     xlab="Observed value",
     ylab="Predicted value", main = "MLR")
abline(a=0, b=1)

#-------------------------MLR test--------------------------------
lm.predicted.y.test <- predict(compounds.lm,newdata = multi.regression.x.test)
lm.predicted.y.test

cor(preprocessed.y.test, lm.predicted.y.test)
lm.r2.test <- cor(preprocessed.y.test, lm.predicted.y.test)**2
lm.r2.test

plot(preprocessed.y.test, lm.predicted.y.test,
     xlab="Observed value",
     ylab="Predicted value", main = "MLR test")
abline(a=0, b=1)


#------------------------PLS training------------------------------------
library(pls)

compounds.plsr <- plsr(preprocessed.y~., data=multi.regression.compounds.train, validation="CV")
summary(compounds.plsr)
plot(compounds.plsr, "validation")

ncomp.onesigma <- selectNcomp(compounds.plsr, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma

predict(compounds.plsr)[, , ncomp.onesigma]
plsr.predicted.y <- predict(compounds.plsr)[, , ncomp.onesigma]
plsr.r2 <- cor(multi.regression.compounds.train[,c(1)], plsr.predicted.y)**2
plsr.r2

plot(multi.regression.compounds.train[,c(1)], plsr.predicted.y,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLSR")
abline(a=0, b=1)

compounds.plsr$coefficients[, , ncomp.onesigma]

#-------------------------pls test--------------------------------
pls.predicted.y.test <- predict(compounds.plsr,newdata = multi.regression.x.test)[,, ncomp.onesigma]
cor(preprocessed.y.test, pls.predicted.y.test)
pls.r2.test <- cor(preprocessed.y.test, pls.predicted.y.test)**2
pls.r2.test

plot(preprocessed.y.test, pls.predicted.y.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS test")
abline(a=0, b=1)

#--------------------------------------PLS-VIP training-------------------------
library(plsVarSel)

vip.selected <- bve_pls(preprocessed.y.train, multi.regression.x.train, ncomp = ncomp.onesigma, VIP.threshold = 0.8) 
vip.selected

x.vip.train <- multi.regression.x.train[,c(vip.selected$bve.selection)]
x.vip.test <- multi.regression.x.test[,c(vip.selected$bve.selection)]
multi.regression.compounds.train.vip <- cbind.data.frame(preprocessed.y.train,x.vip.train)
multi.regression.compounds.test.vip <- cbind.data.frame(preprocessed.y.test,x.vip.test)

compounds.plsr.vip <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.vip, validation="CV")
summary(compounds.plsr.vip)
plot(compounds.plsr.vip, "validation")

ncomp.onesigma.vip <- selectNcomp(compounds.plsr.vip, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.vip

predict(compounds.plsr.vip)[, , ncomp.onesigma.vip]
plsr.predicted.y.vip <- predict(compounds.plsr.vip)[, , ncomp.onesigma.vip]
plsr.r2.vip <- cor(preprocessed.y.train, plsr.predicted.y.vip)**2
plsr.r2.vip

plot(preprocessed.y.train, plsr.predicted.y.vip,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-VIP")
abline(a=0, b=1)

compounds.plsr.vip$coefficients[, , ncomp.onesigma.vip]

#-------------------------pls-vip test--------------------------------
pls.predicted.y.vip.test <- predict(compounds.plsr.vip,newdata = multi.regression.x.test)[,, ncomp.onesigma.vip]
cor(preprocessed.y.test, pls.predicted.y.vip.test)
pls.r2.vip.test <- cor(preprocessed.y.test, pls.predicted.y.vip.test)**2
pls.r2.vip.test

plot(preprocessed.y.test, pls.predicted.y.vip.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-VIP test")
abline(a=0, b=1)


#-----------------------------Uninformative Variable Elimination in PLS training----------------------
mcuve <- mcuve_pls(preprocessed.y.train, multi.regression.x.train, ncomp = ncomp.onesigma, N = 3)
mcuve

x.uve.train <- multi.regression.x.train[,c(mcuve$mcuve.selection)]
x.uve.test <- multi.regression.x.test[,c(mcuve$mcuve.selection)]
multi.regression.compounds.train.uve <- cbind.data.frame(preprocessed.y.train,x.uve.train)
multi.regression.compounds.test.uve <- cbind.data.frame(preprocessed.y.test,x.uve.test)

compounds.plsr.uve <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.uve, validation="CV")
summary(compounds.plsr.uve)
plot(compounds.plsr.uve, "validation")

ncomp.onesigma.uve <- selectNcomp(compounds.plsr.uve, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.uve

predict(compounds.plsr.uve)[, , ncomp.onesigma.uve]
plsr.predicted.y.uve <- predict(compounds.plsr.uve)[, , ncomp.onesigma.uve]
plsr.r2.uve <- cor(preprocessed.y.train, plsr.predicted.y.uve)**2
plsr.r2.uve

plot(preprocessed.y.train, plsr.predicted.y.uve,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-UVE")
abline(a=0, b=1)

compounds.plsr.uve$coefficients[, , ncomp.onesigma.uve]

#-------------------------pls-uve test--------------------------------
pls.predicted.y.uve.test <- predict(compounds.plsr.uve,newdata = multi.regression.x.test)[,, ncomp.onesigma.uve]
cor(preprocessed.y.test, pls.predicted.y.uve.test)
pls.r2.uve.test <- cor(preprocessed.y.test, pls.predicted.y.uve.test)**2
pls.r2.uve.test

plot(preprocessed.y.test, pls.predicted.y.uve.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-uve test")
abline(a=0, b=1)


#-----------------------------GAPLS training-----------------------------------
library(gaselect)

n = ncol(multi.regression.x.train)
l = cbind(ncomp.onesigma + 5,n)

ctrl <- genAlgControl(populationSize = 100, numGenerations = 100, minVariables = 1,
                      maxVariables = min(l), verbosity = 1)

evaluatorRDCV <- evaluatorPLS(numReplications = 2, innerSegments = 5, outerSegments = 3,
                              numThreads = 1)

# Generate demo-data
set.seed(12345)
X <- as.matrix(multi.regression.x.train)
#View(X)
y <- (preprocessed.y.train)

resultRDCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorRDCV, seed = 123)

RD <- subsets(resultRDCV, 1, names = FALSE)
RD$`1`

x.ga.train <- multi.regression.x.train[,c(RD$`1`)]
x.ga.test <- multi.regression.x.test[,c(RD$`1`)]
multi.regression.compounds.train.ga <- cbind.data.frame(preprocessed.y.train,x.ga.train)
multi.regression.compounds.test.ga <- cbind.data.frame(preprocessed.y.test,x.ga.test)

compounds.plsr.ga <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.ga, validation="CV")
summary(compounds.plsr.ga)
plot(compounds.plsr.ga, "validation")

ncomp.onesigma.ga <- selectNcomp(compounds.plsr.ga, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.ga

predict(compounds.plsr.ga)[, , ncomp.onesigma.ga]
plsr.predicted.y.ga <- predict(compounds.plsr.ga)[, , ncomp.onesigma.ga]
plsr.r2.ga <- cor(preprocessed.y.train, plsr.predicted.y.ga)**2
plsr.r2.ga

plot(preprocessed.y.train, plsr.predicted.y.ga,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-GA")
abline(a=0, b=1)

compounds.plsr.ga$coefficients[, , ncomp.onesigma.ga]

#-------------------------pls-GA test--------------------------------
pls.predicted.y.ga.test <- predict(compounds.plsr.ga,newdata = multi.regression.x.test)[,, ncomp.onesigma.ga]
cor(preprocessed.y.test, pls.predicted.y.ga.test)
pls.r2.ga.test <- cor(preprocessed.y.test, pls.predicted.y.ga.test)**2
pls.r2.ga.test

plot(preprocessed.y.test, pls.predicted.y.ga.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-GA test")
abline(a=0, b=1)

#---------------------summary-----------------------
summary.train <-cbind(lm.r2,plsr.r2,plsr.r2.uve,plsr.r2.vip,plsr.r2.ga)
View(summary.train)
summary.test <- cbind(lm.r2.test,pls.r2.test,pls.r2.uve.test,pls.r2.vip.test,pls.r2.ga.test)
View(summary.test)

