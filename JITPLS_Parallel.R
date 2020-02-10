library(readr) # データ読み込み
library(dplyr) # データ操作一般
library(assertr) # データのチェック
library(rsample)

#Packages installation
library(genalg)
library(pls)
library(e1071)
library(kernlab)
library(iml)
library(devtools)

pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#if you want to change the number of threads for the calculation, please change the value "detectCores()"
registerDoParallel(makeCluster(detectCores()))


#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)
#View(compounds)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]
#View(complete.compounds)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]
y

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(2:16)]

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#-----------remove columns with too high cov-----------
library(caret)

df2 <- as.matrix(cor(x))
df2
hc = findCorrelation(df2, cutoff = .6, exact = FALSE) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = x[,-c(hc)]

x <- reduced_Data
x

#---------時系列データ-----------
delay <- 1
range <- 5

stv <- ncol(x)
cd <- matrix(0, ncol = ncol(x)*range, nrow = nrow(x))

cd <- foreach(i = 1:stv, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% {
 foreach(j = 1:range, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% { 
   cd[, c(i + ncol(x)*(j - 1))] <- mutate(x, s = lag(x[,c(i)], delay*(j-1)))[,c(ncol(x) + 1)]
 }
}

colnames(cd) <- paste0("x", 1:ncol(cd))

#------------recollect data with delay---------------------
compounds.rev <- cbind(y, cd)
is.completes.rev <- complete.cases(compounds.rev)
complete.compounds <- compounds.rev[is.completes.rev, ]
View(complete.compounds)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]
y

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(2:66)]
#View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
#View(x)

#-----------standarization of y------------------------
preprocessed.y <- (y - mean(y)) / sd(y)
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x------------------------
apply(x, 2, mean)
apply(x, 2, sd)
preprocessed.x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
#View(preprocessed.x)

#-----------x converted into data frame type for machine learning-----------
#class(preprocessed.x)
preprocessed.x <- data.frame(preprocessed.x)

#-----------compare the number of columns and rows--------------
ncol(preprocessed.x)
nrow(preprocessed.x)


#-----------pick up columns if needed---------------------------
multi.regression.x <- x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(y, x)

#--------------------divide into test and training data----------------------
train_size = 0.9

n = nrow(multi.regression.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[c(1:round(n*train_size)), ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-c(1:round(n*train_size)), ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)
#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

nsamp <- 200

# Euclidean distances
my_distances1 <- dist(multi.regression.x, method = "minkowski")
m <- as.matrix(my_distances1)

query <- matrix(0, nrow = nsamp, ncol = nrow(multi.regression.x.test.s.t))
#View(query)

for(k in 1:nrow(multi.regression.x.test.s.t)){
  query[,c(k)] <- as.matrix(order(m[-c(nrow(multi.regression.x.train.s.t)+1:nrow(multi.regression.x)),c(nrow(multi.regression.x.train.s.t)+k)])[1:nsamp], nrow=1)
}

#View(query)

git <- foreach(j = 1:nrow(multi.regression.compounds.test.s.t), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% {
dat <- matrix(multi.regression.x.test.s.t[c(j),], ncol = ncol(multi.regression.x.test.s.t)) #蝗ｺ螳壹☆繧区擅莉ｶ縺後≠繧句?ｴ蜷医?ｯ縺薙％縺ｧ髯､蜴ｻ

multi.regression.x.train.s.t.git <- multi.regression.x.train.s.t[query[,c(j)],]
x.sds <- apply(multi.regression.x.train.s.t.git, 2, sd)
sd.is.not.0 <- x.sds != 0
multi.regression.x.train.s.t.git <- multi.regression.x.train.s.t.git[, sd.is.not.0]
dat <- matrix(dat[sd.is.not.0], ncol = ncol(multi.regression.x.train.s.t.git))
multi.regression.x.train.s.t.git <- scale(multi.regression.x.train.s.t.git) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(multi.regression.x.train.s.t.git, "scaled:center") 
col_stddevs_train <- attr(multi.regression.x.train.s.t.git, "scaled:scale")
dat <- scale(dat, center = col_means_train, scale = col_stddevs_train)

preprocessed.y.train.git <- preprocessed.y.train[query[,c(j)]]
centy <- mean(preprocessed.y.train.git)
sdy <- sd(preprocessed.y.train.git)
preprocessed.y.train.git <-(preprocessed.y.train.git - centy) / sdy
multi.regression.compounds.train.s.t.git <- as.data.frame(cbind(preprocessed.y.train.git, multi.regression.x.train.s.t.git))

model <- plsr(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, validation = "CV")
ncompopt <- (order(as.matrix(model$validation$PRESS)))[1]
predic <- sdy*predict(model, newdata = dat)[,,ncompopt] + centy
data.frame(preprocessed.y.test[j], predic)
}

plot(git[,c(1)], git[,c(2)])
r2jit <- 1 - sum((git[,c(1)] - git[,c(2)])^2) / sum((mean(preprocessed.y.train) - git[,c(2)])^2)
r2jit             
