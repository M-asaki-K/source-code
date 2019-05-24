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
x.0 <- complete.compounds[,c(11:14)]
x.2 <- apply(x.0,2,function(x.0) {x.0**2}) #add nonlinear columns
x.3 <- apply(x.0,2,function(x.0) {x.0**3}) #add nonlinear columns

x <- cbind.data.frame(x.0)


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
train_size = 1

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

#-------------Installing SVM library---------------------------------
library(e1071)
library(kernlab)
library(iml)
library(devtools)
library(boot)

#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

#------------------------------------outlier sample elimination-----------------------------------------------------------------
zs <- as.data.frame(cbind(preprocessed.y.train,multi.regression.x.train.s.t))

i = 0
for(i in 0:1){

   gam <- matrix(data = 0, nrow = 31, ncol = 1)
   for(k in -20:10){
      rbf <- rbfdot(sigma = 2^k)
      rbf
      
      asmat <- as.matrix(zs[, -c(1)])
      asmat
      
      kern <- kernelMatrix(rbf, asmat)
      sd(kern)
      gam[c(k + 21),] <- sd(kern)
   }
   
   
   hakata <- which.max(gam)
   
   obj.se.t <- tune.svm(zs[,-c(1)], zs[,c(1)], gamma = 2^(hakata - 21), cost = 3, epsilon = 2^(-10:0))
   obj.sc.t <- tune.svm(zs[,-c(1)], zs[,c(1)], gamma = 2^(hakata - 21), cost = 2^(-5:10), epsilon = obj.se.t$best.parameters[,c(3)])
   obj.s.t <- tune.svm(zs[,-c(1)], zs[,c(1)], gamma = 2^(-20:10), cost = obj.sc.t$best.parameters[,c(2)], epsilon = obj.se.t$best.parameters[,c(3)])
   
theta <- function(zs,ind){
   
compounds.svr.s.t <- svm(zs[ind,-c(1)], zs[ind,c(1)],gammma = obj.s.t$best.parameters[,c(1)], cost = obj.s.t$best.parameters[,c(2)], epsilon = obj.s.t$best.parameters[,c(3)])
abs(predict(compounds.svr.s.t, newdata = zs[ind,-c(1)]) - median(predict(compounds.svr.s.t, newdata = zs[ind,-c(1)])))}

theta(zs)

sample <- boot(data = zs,statistic = theta, sim = "ordinary", stype = "i", R =10)
sample
threshold <- 3*1.4826*median(sample$t)
threshold

new.data <- zs[,-c(1)]
compounds.svr.s.t <- svm(zs[,-c(1)], zs[,c(1)],gammma = obj.s.t$best.parameters[,c(1)], cost = obj.s.t$best.parameters[,c(2)], epsilon = obj.s.t$best.parameters[,c(3)])
new.fit <- predict(compounds.svr.s.t, new.data)

nuke.fun <- function(zs, ind, i.pred, fit.pred, x.pred)
{
   lm.b <- svm(zs[ind,-c(1)], zs[ind,c(1)],gammma = obj.s.t$best.parameters[,c(1)], cost = obj.s.t$best.parameters[,c(2)], epsilon = obj.s.t$best.parameters[,c(3)])
   pred.b <- predict(lm.b, x.pred)
   pred.b
}

nuke.boot <- boot(zs, nuke.fun, R = 10, m = 1, 
                  fit.pred = new.fit, x.pred = new.data)

plot(t(apply(nuke.boot$t, 2, median)),zs[,c(1)],
     xlab="Observed value",
     ylab="Predicted value", main = "SVM test")
abline(a=0, b=1)
err <- abs(t(apply(nuke.boot$t, 2, median)) - zs[,c(1)])

if(max(err) < threshold){
   i = 1
}
else{zs <- zs[err < threshold,]
View(zs)
   i = 0
}
}
View(zs) #外れ値を除いたデータを表示する
