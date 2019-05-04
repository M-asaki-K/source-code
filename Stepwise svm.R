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
x.2 <- apply(x.0,2,function(x.0) {x.0**2}) #add nonlinear columns
x.3 <- apply(x.0,2,function(x.0) {x.0**3}) #add nonlinear columns

x <- cbind.data.frame(x.0,x.2,x.3)


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

#-------------Installing SVM library---------------------------------
library(e1071)

#--------------------------------Stepwise variables elimination in SVM---------------------------
#Generating initial dataset(train and test)
multi.regression.compounds.train.s <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s <- multi.regression.x.train[,]
multi.regression.compounds.test.s <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s <- multi.regression.x.test[,]

#Definition of best performance data for each cases
best.performance.cv <- matrix(data = 0, nrow = ncol(multi.regression.x.train), ncol = 2)
best.performance.cv[,c(2)] <- seq(1:ncol(multi.regression.x.train))
best.performance.cv
variables.step <- colnames(preprocessed.x)

#Stepwise variables elimination
for(j in 1:(ncol(preprocessed.x) - 15)){ #the process stops when variable number becomes 1

#SVM hyperparameter tuning  
obj.s <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(-3:-1), cost = 2^(0:3), epsilon = 2^(-3:-1))
obj.s$best.model
compounds.svr.s <- svm(multi.regression.x.train.s,preprocessed.y.train,gammma = obj.s$best.parameters[,c(1)], cost = obj.s$best.parameters[,c(2)], epsilon = obj.s$best.parameters[,c(3)])
summary(compounds.svr.s)
best.performance.cv[c(j),c(1)] <- obj.s$best.performance #input the best performance for each steps

variable.importance.sum <- matrix(data = 0, nrow = ncol(multi.regression.x.train.s), ncol = 1)
rownames(variable.importance.sum) <- colnames(multi.regression.x.train.s)

#calculating variable importance for each factor
for (i in 1:ncol(multi.regression.x.train.s)){
  pdp.data <- matrix(data = 0, nrow = 10, ncol = ncol(multi.regression.x.train.s))
  d <- seq(-1,1,length = 10)
  d
  pdp.data[,c(i)]<- d
  variable.importance <- max(predict(compounds.svr.s,newdata = pdp.data)) - min (predict(compounds.svr.s,newdata = pdp.data))
  pdp.data <- matrix(data = 0, nrow = 10, ncol = ncol(multi.regression.x.train.s))
  variable.importance.sum[c(i),] <- variable.importance
  
}
#eliminate variable with lowest importance
eliminate <- which.min(variable.importance.sum)

multi.regression.compounds.train.s <- cbind(preprocessed.y.train, multi.regression.x.train.s[,-c(eliminate)])
multi.regression.x.train.s <- multi.regression.x.train.s[,-c(eliminate)]
multi.regression.compounds.test.s <- cbind(preprocessed.y.test, multi.regression.x.test.s[,-c(eliminate)])
multi.regression.x.test.s <- multi.regression.x.test.s[,-c(eliminate)]

#output selected variables for each steps
variables.step <- cbind(variables.step, colnames(multi.regression.x.train.s))
View(variables.step)
View(best.performance.cv)
}

#output the selected variables
colnames(multi.regression.x.train.s)

#generating SVM model with selected variables
obj.s.t.t <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(-3:-1), cost = 2^(0:4), epsilon = 2^(-3:-1))
obj.s.t.t$best.model
compounds.svr.s.t.t <- svm(multi.regression.x.train.s,preprocessed.y.train,gammma = obj.s.t.t$best.parameters[,c(1)], cost = obj.s.t.t$best.parameters[,c(2)], epsilon = obj.s.t.t$best.parameters[,c(3)])
summary(compounds.svr.s.t.t)
obj.s.t.t$best.performance

#testing the model accuracy
svm.predicted.y.test.s.t.t <- predict(compounds.svr.s.t.t, newdata = multi.regression.x.test.s)
plot(preprocessed.y.test, svm.predicted.y.test.s.t.t,
     xlab="Observed value",
     ylab="Predicted value", main = "SVM test")
abline(a=0, b=1)

svm.r2.test.s.t.t <- cor(preprocessed.y.test,svm.predicted.y.test.s.t.t)**2
svm.r2.test.s.t.t

#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

#model optimization with all variables
obj.s.t <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s.t, gamma = 2^(-3:-1), cost = 2^(0:3), epsilon = 2^(-3:-1))
compounds.svr.s.t <- svm(multi.regression.x.train.s.t,preprocessed.y.train,gammma = obj.s.t$best.parameters[,c(1)], cost = obj.s.t$best.parameters[,c(2)], epsilon = obj.s.t$best.parameters[,c(3)])

variable.importance.sum.t <- matrix(data = 0, nrow = ncol(multi.regression.x.train.s.t), ncol = 1)
rownames(variable.importance.sum.t) <- colnames(multi.regression.x.train.s.t)

#calculating variable importance
for (i in 1:ncol(multi.regression.x.train.s.t)){
  pdp.data.t <- matrix(data = 0, nrow = 10, ncol = ncol(multi.regression.x.train.s.t))
  d <- seq(-1,1,length = 10)
  d
  pdp.data.t[,c(i)]<-d
  variable.importance.t <- max(predict(compounds.svr.s.t,newdata = pdp.data.t)) - min (predict(compounds.svr.s.t,newdata = pdp.data.t))
  pdp.data.t <- matrix(data = 0, nrow = 10, ncol = ncol(multi.regression.x.train.s.t))
  variable.importance.sum.t[c(i),] <- variable.importance.t
  
}

#elimination of variables with low importance
for(k in 1:(ncol(preprocessed.x) - 15)){ 
  eliminate.t <- which.min(variable.importance.sum.t)
  multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train.s.t[,-c(eliminate.t)])
  multi.regression.x.train.s.t <- multi.regression.x.train.s.t[,-c(eliminate.t)]
  multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test, multi.regression.x.test.s.t[,-c(eliminate.t)])
  multi.regression.x.test.s.t <- multi.regression.x.test.s.t[,-c(eliminate.t)]
  
  variable.importance.sum.t <- variable.importance.sum.t[-c(s),]
  variable.importance.sum.t <- as.matrix(variable.importance.sum.t)
}
variable.importance.sum.t
View(multi.regression.x.train.s.t)

#generating SVM model with selected variables
obj.s.t.e <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s.t, gamma = 2^(-3:-1), cost = 2^(0:4), epsilon = 2^(-3:-1))
obj.s.t.e$best.model
compounds.svr.s.t.e <- svm(multi.regression.x.train.s.t,preprocessed.y.train,gammma = obj.s.t.e$best.parameters[,c(1)], cost = obj.s.t.e$best.parameters[,c(2)], epsilon = obj.s.t.e$best.parameters[,c(3)])
summary(compounds.svr.s.t.e)
obj.s.t.e$best.performance

#testing the model accuracy
svm.predicted.y.test.s.t.e <- predict(compounds.svr.s.t.e, newdata = multi.regression.x.test.s.t)
plot(preprocessed.y.test, svm.predicted.y.test.s.t.e,
     xlab="Observed value",
     ylab="Predicted value", main = "SVM test")
abline(a=0, b=1)

svm.r2.test.s.t.e <- cor(preprocessed.y.test,svm.predicted.y.test.s.t.e)**2
svm.r2.test.s.t.e

