library("keras")

#-----------read csv file as compounds--------------
#path <- file.choose()
compounds <- read.csv("../input/traintakeda/train.csv")
#compounds <- read.csv(path)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,-c(1)]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]

#-----------read csv file as compounds--------------
#path.t <- file.choose()
compounds.t <- read.csv("../input/testtakeda/test.csv")
#compounds.t <- read.csv(path.t)

#-----------remove some columns if needed--------------
trimed.compounds.t <- compounds.t[,-c(1)]

#-----------select rows without empty cells---------
is.completes.t <- complete.cases(trimed.compounds.t)
is.completes.t

complete.compounds.t <- trimed.compounds.t[is.completes.t,]

#-----------select x from the dataset-----------------
x.0 <- complete.compounds[,-c(1)]
x.0.t <- complete.compounds.t[,]
x.sds <- apply(x.0, 2, sd)
x.sds.t <- apply(x.0.t, 2, sd)

sd.is.not.0 <- x.sds != 0 
sd.is.0 <- x.sds == 0
sd.is.not.0.t <- x.sds.t != 0
sd.is.0.t <- x.sds.t == 0

x.0.eda <- x.0[, c(sd.is.0.t)]
x.0.eda.t <- x.0.t[, c(sd.is.0)]

row_sub = apply(x.0.eda[,-c(16)], 1, function(row) all(row ==0 ))
row_sub_2 = x.0.eda[,c(16)] ==1
row_sub.t = apply(x.0.eda.t, 1, function(row) all(row ==0 ))
row_sub.t[row_sub.t == FALSE]

sd.is.not.0.f <- sd.is.not.0.t&sd.is.not.0
row_sub.f <- row_sub&row_sub_2
row_sub.t.f <- row_sub.t
x.0 <- x.0[row_sub.f, (sd.is.not.0.f)]
x.0.t <- x.0.t[, (sd.is.not.0.f)]

x.sds <- apply(x.0, 2, sd)
x.sds.t <- apply(x.0.t, 2, sd)

sd.is.not.0 <- x.sds != 0 
sd.is.not.0.t <- x.sds.t != 0
sd.is.not.0.f <- sd.is.not.0.t&sd.is.not.0
x.0 <- x.0[, (sd.is.not.0.f)]
x.0.t <- x.0.t[, (sd.is.not.0.f)]

nrow(x.0)
ncol(x.0)
nrow(x.0.t)
ncol(x.0.t)

x.1 <- rbind(x.0,x.0.t)
pc = 500
rpca = prcomp(x = x.1,scale=F)
x.0.s <- as.data.frame(rpca$x[c(1:nrow(x.0)), c(1:pc)])
x.0.t.s <- as.data.frame(rpca$x[-c(1:nrow(x.0)), c(1:pc)])

#----------------------mahalanobis distance calculation--------------------------------------
#Mahalanobis <- function(dat,                 # 基準群のデータ行列
#                        x)                      # 所属確率を計算するデータ行列
#{
#  dat <- subset(dat, complete.cases(dat))      # 欠損値を持つケースを除く
#  n <- nrow(dat)                               # ケース数
#  p <- ncol(dat)                               # 変数の個数
#  ss <- var(dat)*(n-1)/n                       # 分散・共分散行列
#  inv.ss <- solve(ss)                  # 分散共分散行列の逆行列
#  m <- colMeans(dat)                   # 各変数の平均値
#  dif <- t(t(x)-m)                     # 平均値からの偏差
#  d2 <- apply(dif, 1,
#              function(z) z %*% inv.ss %*% z) # マハラノビスの平方距離
#  P <- pchisq(d2, p, lower.tail=FALSE) # 所属確率
#  return(data.frame(d2=d2, P=P))
#}

#model.threshold <- max(Mahalanobis(x.0.s, x.0.s))
#model.threshold

#--------------------divide into test and training data----------------------
train_size = 0.6

n = nrow(x.0.s)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

train_data <- cbind.data.frame(x.0.s[perm, ])
test_data <- cbind.data.frame(x.0.s[-perm, ])
nrow(test_data)
#-----------select y from the dataset------------------
train_labels <- complete.compounds[row_sub.f, c(1)]
train_labels <- train_labels[perm]
test_labels <- complete.compounds[row_sub.f, c(1)]
test_labels <- test_labels[-perm]
test_labels

# Test data is *not* used when calculating the mean and std.

# Normalize training data
train_data <- scale(train_data) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(train_data, "scaled:center") 
col_stddevs_train <- attr(train_data, "scaled:scale")
test_data <- scale(test_data, center = col_means_train, scale = col_stddevs_train)

train_data[1, ] # First training sample, normalized

#-----------select x from the dataset-----------------
test_data.t <- cbind.data.frame(x.0.t.s[,])
#model.threshold.t <- Mahalanobis(x.0.s, test_data.t)[, c(1)] < model.threshold
#model.threshold.t

# Test data is *not* used when calculating the mean and std.
# Use means and standard deviations from training set to normalize test set
test_data.t <- scale(test_data.t, center = col_means_train, scale = col_stddevs_train)

build_model <- function() {
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 1000, activation = "relu",
                input_shape = dim(train_data)[2]) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.8) %>%
    layer_dense(units = 300, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.8) %>%
#    layer_dense(units = 200, activation = "relu") %>%
#    layer_batch_normalization() %>%
#    layer_dropout(rate = 0.8) %>%    
#    layer_dense(units = 100, activation = "relu") %>%
#    layer_batch_normalization() %>%
#    layer_dropout(rate = 0.75) %>%
    layer_dense(units = 1)
  
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adadelta(lr = initial_lr),
    metrics = list("mean_absolute_error")
  )
  
  model
}

epochs <- 2000
batch_size <- 128
initial_lr<-0.2
decay<-2
period<-25


model <- build_model()
model %>% summary()

# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)    

lr_schedule<-function(epoch,lr) {
  lr=initial_lr/decay^((epoch-1)%%period)+1e-5
  lr
}
cb_lr<-callback_learning_rate_scheduler(lr_schedule)
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 100)
#rr <- callback_learning_rate_scheduler(function(epochs){0.01 / epochs})

# Fit the model and store training stats
history <- model %>% fit(
  train_data,
  train_labels,
  epochs = epochs,
  validation_split = 0.2,
  verbose = 0,
  batch_size = batch_size,
  callbacks = list(print_dot_callback, cb_lr, early_stop)
)

plot(history)

train_predictions <- model %>% predict(train_data)
plot(train_predictions, train_labels)
r2train <- cor(train_predictions, train_labels)**2
r2train

test_predictions <- model %>% predict(test_data)
plot(test_predictions, test_labels)
r2test <- cor(test_predictions, test_labels)**2
r2test

test_predictions.t <- model %>% predict(test_data.t)
test_predictions.t <- as.matrix(test_predictions.t)
#test_predictions.t[row_sub.t == FALSE,] <- mean(complete.compounds[,c(1)]) found to be not effective for this case
#test_predictions.t[-(model.threshold.t)] <- mean(complete.compounds[,c(1)])
write.csv(test_predictions.t[1:13732])
