# another realistic example: wavelenght selection for PLS on NIR data
## Not run:
library(genalg)
library(pls)
library(e1071)
library(kernlab)
library(iml)
library(devtools)

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

#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

Xtrain <- multi.regression.x.train.s.t
ytrain <- preprocessed.y.train
datasum <- as.data.frame(cbind(ytrain,Xtrain))

numberOfWavelenghts = ncol(data0) - 1

evaluateNIR <- function(chromosome=c()) {
  returnVal = 100
  minLV = 2
  if (sum(chromosome) < minLV) {
    returnVal
  } else {
    xtrain = Xtrain[,chromosome == 1];
    datasum = as.data.frame(cbind(ytrain, xtrain))
    
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
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
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
    parms
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:max(df2$fold), .combine = rbind, .inorder = FALSE) %dopar% {
        deve <- df2[df2$fold != j, ]
        test <- df2[df2$fold == j, ]
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
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
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
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
      }
}

monitor <- function(obj) {
  minEval = min(obj$evaluations);
  filter = obj$evaluations == minEval;
  bestObjectCount = sum(rep(1, obj$popSize)[filter]);
  # ok, deal with the situation that more than one object is best
  if (bestObjectCount > 1) {
    bestSolution = obj$population[filter,][1,];
  } else {
    bestSolution = obj$population[filter,];
  }
  outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
                     " Best (Error=", minEval, "): ", sep="");
  for (var in 1:length(bestSolution)) {
    outputBest = paste(outputBest,
                       bestSolution[var], " ",
                       sep="");
  }
  outputBest = paste(outputBest, "\n", sep="");
  cat(outputBest);
}

nir.results = rbga.bin(size=numberOfWavelenghts, zeroToOneRatio=10,
                       evalFunc=evaluateNIR, monitorFunc=monitor,
                       popSize=20, iters= 20, verbose=TRUE)

nir.results
