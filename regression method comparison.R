#-----------------dataset--------------------
b <- c(0,0.2,0.4,0.6,0.8,1)
d <- c(0.1,0.2,0.4,0.6,0.7,0.8)
b <- as.data.frame(b)
d <- as.data.frame(d)

h <- as.data.frame(c(d,b))

#-----------linear regression----------------
linear <- lm(d~.,h)
summary(linear)
pred.linear <- predict.lm(linear,b)
plot(h,pred.linear)
View(pred.linear)

#------------gaussian----------------
library(GPfit)
dim(b)
dim(d)
x = h[,c(2)];
y = h[,c(1)];
GPmodel = GP_fit(x,y);
print(GPmodel)

m <- c(-1,-0.5,-0.2,0.2,0.4,0.6,0.8,1.2,1.6,2) # valiable
M.data <-cbind.data.frame(m) 
M.data

GPprediction = predict.GP(GPmodel,xnew=M.data);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
plot(m,yhat)

sigma <- 3*(mse**0.5) # 3 sigma
f <- cbind(yhat, m, sigma)
f <- as.data.frame(f)

library(reshape2)
library(ggplot2)
library(ggsci)

g <- ggplot(f, aes(x = m, y = yhat))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin = yhat - sigma, ymax = yhat + sigma, width = 0.03))
g <- g + scale_fill_nejm()
g <- g + xlim(-1,2) + ylim(-0.5,1.5)

k <- cbind(m,yhat,sigma)
View(k)
plot(g)

#-------------------randomforest-------------------
b <- c(0,0.2,0.4,0.6,0.8,1,-1,-0.5,0,0.2,0.4,0.6,0.8,1,1.5,2)
d <- c(0.1,0.2,0.4,0.6,0.7,0.8,-0.7,-0.2,0.1,0.2,0.4,0.6,0.7,0.8,1.2,1.6)
b <- as.data.frame(b)
d <- as.data.frame(d)

h <- as.data.frame(c(d,b))
set.seed(131)
rand <- randomForest(d ~ ., data=h[-c(7:16)], mtry=1,
                         importance=TRUE, na.action=na.omit)
print(rand)

pred.rand <- predict(rand,h[c(7:16),])
l <- cbind(h[c(7:16),c(2)],pred.rand)
View(l)

#----------------------mahalanobis distance calculation--------------------------------------
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

#使用例

dat <- h
b1 <- c(0,0.2,0.4,0.6,0.8,1,-1,-0.5,0,0.2,0.4,0.6,0.8,1,1.5,2)
d1 <- c(0.1,0.2,0.4,0.6,0.7,0.8,-0.7,-0.2,0.1,0.2,0.4,0.6,0.7,0.8,1.2,1.6)
b1 <- as.data.frame(b)
d1 <- as.data.frame(d)

h1 <- as.data.frame(c(d1,b1))

Mahalanobis(dat, h1)	

