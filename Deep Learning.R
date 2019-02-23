# Installation
cran <- getOption("repos")
cran["dmlc"] <- "https://apache-mxnet.s3-accelerate.dualstack.amazonaws.com/R/CRAN/"
options(repos = cran)
install.packages("mxnet")

library(mxnet)
a <- mx.nd.ones(c(2,3), ctx = mx.cpu())
b <- a * 2 + 1
b

scaled.compounds <- cbind(nominal.y, preprocessed.x)
#View(scaled.compounds)

library(mxnet)

train_size = 0.2

n = nrow(scaled.compounds)
# 1 ～ n から無作為に n * train_size 個を抽出
perm = sample(n, size = round(n * train_size))

# 学習用データ
train <- scaled.compounds[perm, ]
# 評価用データ
test <-scaled.compounds[-perm, ]
View(train)

# 学習用入力データ
train.x <- data.matrix(train[,-c(1)])
# 学習用ラベルデータ（0 ～ 2）
train.y <- as.numeric(train[,c(1)])-1
train.y

# 評価用入力データ
test.x <- data.matrix(test[,-c(1)])
# 評価用ラベルデータ（1 ～ 3）
test.y <- as.numeric(test[,c(1)])
test.y

mx.set.seed(0)

# 学習
model.cnnl <- mx.mlp(train.x, train.y, 
                hidden_node = c(300,300,300),
                out_node = 2,dropout = 0.6,
                num.round = 100,
                learning.rate = 0.06,
                array.batch.size = 10,
                activation = 'relu',
                array.layout = 'rowmajor',
                eval.metric = mx.metric.accuracy)

# 評価
pred <- predict(model.cnnl, test.x, array.layout = 'rowmajor')

# 評価用データの分類結果（1 ～ 3）
pred.y <- max.col(t(pred))
pred.y

# 評価データの正解率を算出
acc <- sum(pred.y == test.y) / length(pred.y)

print(acc)

table(pred.y,test.y)

#異常値の正解率を算出
n <- as.matrix(table(pred.y,test.y))
n <- n[2,2]
accn <- sum(n) / sum(test.y == 2)
print(accn)

help(mxnet)
