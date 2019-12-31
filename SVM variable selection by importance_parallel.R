pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
registerDoParallel(makeCluster(detectCores()))

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
x.0 <- complete.compounds[,-c(1,2)]
x <- cbind.data.frame(x.0)


#-----------remove columns of 0 distribution from x----
x.sds <- apply(x, 2, sd)
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#-----------select y from the dataset------------------
y <- complete.compounds[,c(2)]

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
train_size = 0.9

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

#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

xtrain = multi.regression.x.train.s.t
datasum = multi.regression.compounds.train.s.t

df2 <- datasum
### SPLIT DATA INTO K FOLDS ###
set.seed(2016)
df2$fold <- caret::createFolds(1:nrow(df2), k = 5, list = FALSE)
### PARAMETER LIST ###
cost <- 3
epsilon <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 0)

gam <- matrix(data = 0, nrow = 31, ncol = 1)
for(k in -20:10){
  rbf <- rbfdot(sigma = 2^k)
  rbf
  
  asmat <- as.matrix(xtrain)
  asmat
  
  kern <- kernelMatrix(rbf, asmat)
  sd(kern)
  gam[c(k + 21),] <- sd(kern)
}

hakata <- which.max(gam)

gamma <- hakata - 21
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- df2[df2$fold != j, ]
    test <- df2[df2$fold == j, ]
    mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
  data.frame(parms[i, ], roc)
}

epsilon <- min(result[result[, c(4)] <= (min(result[,c(4)])), c(1)])
cost <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- df2[df2$fold != j, ]
    test <- df2[df2$fold == j, ]
    mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
  data.frame(parms[i, ], roc)
}

cost <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(2)])
gamma <- c(-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
    deve <- df2[df2$fold != j, ]
    test <- df2[df2$fold == j, ]
    mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
  data.frame(parms[i, ], roc)
}

gamma <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(3)])
bestperformance <- min(result[, c(4)])
returnVal = bestperformance
compounds.svr.s.t <- svm(multi.regression.x.train.s.t, preprocessed.y.train, gammma = 2^gamma, cost = 2^cost, epsilon = 2^epsilon)

#------------------------------feature importance calculation----------------------------------------------
mod = Predictor$new(compounds.svr.s.t, data = multi.regression.x.train.s.t, y = preprocessed.y.train)
imp = FeatureImp$new(mod, loss = "mse", compare = "ratio", n.repetitions = 5, parallel = TRUE)
imp$results
plot(imp)

#--------------------variable selection threshold by onesigma method------------------
threshold <- min(imp$results[, c(3)]) + mean(imp$results[, c(4)] - imp$results[, c(2)]) / 3.92*4
acq <- cbind(imp$results[, c(1)], imp$results[, c(3)])
kamakura <- acq[acq[, c(2)] > threshold, c(1)]
kamakura

#------------elimination of variables with low importance------------------------------
  multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train.s.t[, c(kamakura)])
  multi.regression.x.train.s.t <- multi.regression.x.train.s.t[, c(kamakura)]
  multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test, multi.regression.x.test.s.t[, c(kamakura)])
  multi.regression.x.test.s.t <- multi.regression.x.test.s.t[, c(kamakura)]

#---------------------generating SVM model with selected variables----------------------
  xtrain = multi.regression.x.train.s.t
  datasum = multi.regression.compounds.train.s.t
  
  df2 <- datasum
  ### SPLIT DATA INTO K FOLDS ###
  set.seed(2016)
  df2$fold <- caret::createFolds(1:nrow(df2), k = 5, list = FALSE)
  ### PARAMETER LIST ###
  cost <- 3
  epsilon <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 0)
  
  gam <- matrix(data = 0, nrow = 31, ncol = 1)
  for(k in -20:10){
    rbf <- rbfdot(sigma = 2^k)
    rbf
    
    asmat <- as.matrix(Xtrain)
    asmat
    
    kern <- kernelMatrix(rbf, asmat)
    sd(kern)
    gam[c(k + 21),] <- sd(kern)
  }
  
  hakata <- which.max(gam)
  
  gamma <- hakata - 21
  parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
  ### LOOP THROUGH PARAMETER VALUES ###
  result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
    c <- parms[i, ]$cost
    g <- parms[i, ]$gamma
    e <- parms[i, ]$epsilon
    ### K-FOLD VALIDATION ###
    out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
      deve <- df2[df2$fold != j, ]
      test <- df2[df2$fold == j, ]
      mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
      pred <- predict(mdl, test)
      data.frame(test[, c(1)], pred)
    }
    ### CALCULATE SVM PERFORMANCE ###
    roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
    data.frame(parms[i, ], roc)
  }
  
  epsilon <- min(result[result[, c(4)] <= (min(result[,c(4)])), c(1)])
  cost <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
  parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
  ### LOOP THROUGH PARAMETER VALUES ###
  result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
    c <- parms[i, ]$cost
    g <- parms[i, ]$gamma
    e <- parms[i, ]$epsilon
    ### K-FOLD VALIDATION ###
    out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
      deve <- df2[df2$fold != j, ]
      test <- df2[df2$fold == j, ]
      mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
      pred <- predict(mdl, test)
      data.frame(test[, c(1)], pred)
    }
    ### CALCULATE SVM PERFORMANCE ###
    roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
    data.frame(parms[i, ], roc)
  }
  
  cost <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(2)])
  gamma <- c(-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
  parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
  ### LOOP THROUGH PARAMETER VALUES ###
  result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
    c <- parms[i, ]$cost
    g <- parms[i, ]$gamma
    e <- parms[i, ]$epsilon
    ### K-FOLD VALIDATION ###
    out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
      deve <- df2[df2$fold != j, ]
      test <- df2[df2$fold == j, ]
      mdl <- e1071::svm(preprocessed.y.train~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
      pred <- predict(mdl, test)
      data.frame(test[, c(1)], pred)
    }
    ### CALCULATE SVM PERFORMANCE ###
    roc <- sum((out[, c(1)] - out[, c(2)])^2) / ncol(datasum)
    data.frame(parms[i, ], roc)
  }
  
  gamma <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(3)])
  bestperformance <- min(result[, c(4)])
  returnVal = bestperformance
  compounds.svr.s.t.e <- svm(multi.regression.x.train.s.t, preprocessed.y.train, gammma = 2^gamma, cost = 2^cost, epsilon = 2^epsilon)

#--------------------------testing the model accuracy------------------------------------------
svm.predicted.y.test.s.t.e <- predict(compounds.svr.s.t.e, newdata = multi.regression.x.test.s.t)
plot(preprocessed.y.test, svm.predicted.y.test.s.t.e,
     xlab="Observed value",
     ylab="Predicted value", main = "SVM test")
abline(a=0, b=1)

svm.r2.test.s.t.e <- cor(preprocessed.y.test,svm.predicted.y.test.s.t.e)**2
svm.r2.test.s.t.e

