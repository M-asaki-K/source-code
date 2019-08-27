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
View(complete.compounds[c(1:5),])

#-----------select x from the dataset-----------------
x <- complete.compounds[,-c(1)]
x <- cbind.data.frame(x)

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
train_size = 0.5

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

#---------------------------stepwise variables selection--------------------------------------------
#Generating initial dataset(train and test)
multi.regression.compounds.train.s <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s <- multi.regression.x.train[,]
multi.regression.compounds.test.s <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s <- multi.regression.x.test[,]

#Definition of best performance data for each cases
best.performance.cv <- matrix(data = 0, nrow = ncol(multi.regression.x.train), ncol = 3)
best.performance.cv[,c(3)] <- seq(1:ncol(multi.regression.x.train))
variables.step <- colnames(preprocessed.x)
importances <- matrix(data = 0, nrow = ncol(multi.regression.x.train), ncol = ncol(multi.regression.x.train))
rownames(importances) <- colnames(preprocessed.x)

#Stepwise variables elimination
for(j in 1:(ncol(preprocessed.x) - 7)){ #the process stops when variable number becomes 6 + 1 (determined from best performance curve)
  
  #SVM hyperparameter tuning  
  #determining initial gamma by maximizing kernel matrix variance
  gam <- matrix(data = 0, nrow = 31, ncol = 1)
  for(k in -20:10){
    rbf <- rbfdot(sigma = 2^k)
    rbf
    
    asmat <- as.matrix(multi.regression.x.train.s)
    asmat
    
    kern <- kernelMatrix(rbf, asmat)
    sd(kern)
    gam[c(k + 21),] <- sd(kern)
  }
  
  hakata <- which.max(gam)
  
  obj.se <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(hakata - 21), cost = 3, epsilon = 2^(-10:0))
  obj.sc <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(hakata - 21), cost = 2^(-5:10), epsilon = obj.se$best.parameters[,c(3)])
  obj.s <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(-20:10), cost = obj.sc$best.parameters[,c(2)], epsilon = obj.se$best.parameters[,c(3)])
  compounds.svr.s <- svm(multi.regression.x.train.s,preprocessed.y.train,gammma = obj.s$best.parameters[,c(1)], cost = obj.s$best.parameters[,c(2)], epsilon = obj.s$best.parameters[,c(3)])
  summary(compounds.svr.s)
  best.performance.cv[c(j), c(1)] <- obj.s$performances[which.min(obj.s$performances[,c(4)]), c(4)]
  best.performance.cv[c(j), c(2)] <- obj.s$performances[which.min(obj.s$performances[,c(4)]), c(5)]
  mod = Predictor$new(compounds.svr.s, data = multi.regression.x.train.s, y = preprocessed.y.train)
  imp = FeatureImp$new(mod, loss = "mse", compare = "ratio", n.repetitions = 10)
  importances[imp$results[,c(1)],c(j)] <- imp$results[,c(3)]
  eliminate <- imp$results[-c(which.min(imp$results[,c(3)])), c(1)]
  eliminate
  
  multi.regression.compounds.train.s <- cbind(preprocessed.y.train, multi.regression.x.train.s[,c(eliminate)])
  multi.regression.x.train.s <- multi.regression.x.train.s[,c(eliminate)]
  multi.regression.compounds.test.s <- cbind(preprocessed.y.test, multi.regression.x.test.s[,c(eliminate)])
  multi.regression.x.test.s <- multi.regression.x.test.s[,c(eliminate)]
  
  #output selected variables for each steps
  variables.step <- cbind(variables.step, colnames(multi.regression.x.train.s))
  View(variables.step)
  View(best.performance.cv)
  View(importances)
}

#imp$results[,c(1)]

#output the selected variables
colnames(multi.regression.x.train.s)

multi.regression.compounds.train.s <- cbind(preprocessed.y.train, multi.regression.x.train[,c(variables.step[c(1:8), c(12)])])
multi.regression.x.train.s <- multi.regression.x.train[,c(variables.step[c(1:8), c(12)])]
multi.regression.compounds.test.s <- cbind(preprocessed.y.test, multi.regression.x.test[,c(variables.step[c(1:8), c(12)])])
multi.regression.x.test.s <- multi.regression.x.test[,c(variables.step[c(1:8), c(12)])]

#generating SVM model with selected variables
gam <- matrix(data = 0, nrow = 31, ncol = 1)
for(k in -20:10){
  rbf <- rbfdot(sigma = 2^k)
  rbf
  
  asmat <- as.matrix(multi.regression.x.train.s.t)
  asmat
  
  kern <- kernelMatrix(rbf, asmat)
  sd(kern)
  gam[c(k + 21),] <- sd(kern)
}

hakata <- which.max(gam)

obj.se.t.t <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(hakata - 21), cost = 3, epsilon = 2^(-10:0))
obj.sc.t.t <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(hakata - 21), cost = 2^(-5:10), epsilon = obj.se.t.t$best.parameters[,c(3)])
obj.s.t.t <- tune.svm(preprocessed.y.train~., data = multi.regression.compounds.train.s, gamma = 2^(-20:10), cost = obj.sc.t.t$best.parameters[,c(2)], epsilon = obj.se.t.t$best.parameters[,c(3)])
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

#--------------------------Bayesian optimization for hyperparameter tuning---------------------------
library(rBayesianOptimization)

Mahalanobis <- function(dat,                 # 基準群のデータ行列
                        x)                      # 所属確率を計算するデータ行列
{
  dat <- subset(dat, complete.cases(dat))      # 欠損値を持つケースを除く
  n <- nrow(dat)                               # ケース数
  p <- ncol(dat)                               # 変数の個数
  ss <- var(dat)*(n-1)/n                       # 分散・共分散行列
  inv.ss <- solve(ss)                  # 分散共分散行列の逆行列
  m <- colMeans(dat)                   # 各変数の平均値
  dif <- t(t(x)-m)                     # 平均値からの偏差
  d2 <- apply(dif, 1,
              function(z) z %*% inv.ss %*% z) # マハラノビスの平方距離
  P <- pchisq(d2, p, lower.tail=FALSE) # 所属確率
  return(data.frame(d2=d2, P=P))
}

dat <- multi.regression.x.train.s[,] #固定する条件がある場合はここで除去
outlier <- max(Mahalanobis(dat, data.frame(dat))[,c(1)])　#maharanobis distによる外れ閾値

#ここからクラスタリングによる外れ値検出
hist <- hclust(dist(multi.regression.x.train.s),method = "ward.D2")
plot(hist)

# 階層型クラスタリングからクラスタ数を手動指定して分類する（コサイン類似度、ウォード法）。
clusnum <- 3
clusters <- cutree(hist, k = clusnum)

# クラスタiの中心点を算出
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

# 全クラスタに適用
centers <- sapply(unique(clusters), clust.centroid, multi.regression.x.train.s, clusters)

# k-meansの実行。sapply()の結果だと行と列が逆なので、転置して引数に与える
km <- kmeans(multi.regression.x.train.s, centers=t(centers)) 
hem <- sqrt(rowSums(multi.regression.x.train.s - fitted(km)) ^ 2)
outlierkm <- mean(hem) + 3*sd(hem)　#kmeansによる外れ値検出

distances <- matrix(data = 0, nrow = clusnum, ncol = 1)
minimumdistance <- function(v){
  for(pp in 1:clusnum){
  distances[c(pp), ] <- sqrt(rowSums(v - km$centers[c(pp),]) ^ 2)
  }
  distances[c(which.min(distances)), ]}

#row 201-203 should be neglected if you do not care about Mahalanobis distance
Gauss_holdout <- function(a, b, c, d, e, f, g){
  daisen <- data.frame((a - 0.5)*6, (b - 0.5)*6, (c - 0.5)*6, (d - 0.5)*6, (e - 0.5)*6, (f - 0.5)*6, (g - 0.5)*6)
{  if(Mahalanobis(dat, data.frame(daisen[,]))[,c(1)] > outlier){ #固定する条件がある場合はここでdaisenから除去
    
    model <- compounds.svr.s.t.t
    Pred <- min(cbind(0, predict(model, data.frame(daisen)))) #最小化問題なら予測値に-1をかける
  list(Score=Pred, Pred=Pred)}
  else{
    if(minimumdistance(daisen) > outlierkm){
      model <- compounds.svr.s.t.t
      Pred <- min(cbind(0, predict(model, data.frame(daisen)))) #最小化問題なら予測値に-1をかける
      list(Score=Pred, Pred=Pred)
    }
    else{
  model <- compounds.svr.s.t.t
  Pred <- predict(model, data.frame(daisen))　#最小化問題なら予測値に-1をかける
  
    list(Score=Pred, Pred=Pred)}}
    }}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(a=c(0,1),b=c(0,1), c=c(0,1), d=c(0,1), e=c(0,1), f=c(0,1), g=c(0,1)), init_points=50, n_iter=3, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)
opt_svm$History

ao <- (opt_svm$Best_Par[1]-0.5)*6
bo <- (opt_svm$Best_Par[2]-0.5)*6
co <- (opt_svm$Best_Par[3]-0.5)*6
do <- (opt_svm$Best_Par[4]-0.5)*6

opt <- cbind(0,ao,0,bo,co,0,do)

#Re-scale of optimized x
A <- as.matrix(t(apply(x[,c(variables.step[c(1:7), c(13)])], 2, sd)))
B <- as.matrix(t(apply(x[,c(variables.step[c(1:7), c(13)])], 2, mean)))
xrev <- t(apply(opt, 1, function(opt){(opt)*A + B}))
colnames(xrev) <- colnames(multi.regression.x.train.s)
View(xrev)

#confimarion of optimized x within the dataset boundaries
xrevs <- opt
#colnames(xrevs) <- kamakura

bmax <- apply(multi.regression.compounds[,c(kamakura)],2,max)
bmax - xrevs > 0

bmin <- apply(multi.regression.compounds[,c(kamakura)],2,min)
xrevs - bmin > 0

#optimized conditions considering mahalanobis distance

path2 <- file.choose()
shiga <- read.csv(path2)
opt2 <- shiga[,-c(1:2,14)]

#Re-scale of optimized x
A <- as.matrix(t(apply(x[,c(kamakura)], 2, sd)))
B <- as.matrix(t(apply(x[,c(kamakura)], 2, mean)))
xrev2 <- t(apply(opt2, 1, function(opt2){(opt2)*A + B}))
colnames(xrev2) <- kamakura
View(xrev2)

min(preprocessed.y)

path <- file.choose()

comp1 <- read.csv(path)

library(rgl)

x = comp1[,c(1)]
y = comp1[,c(2)]
z = comp1[,c(3)]

c = comp1[,c(4)]*(-1)
c = cut(c, breaks=20)
cols = gray.colors(20)[as.numeric(c)]

plot3d(x,y,z,size = 10, col=cols)
