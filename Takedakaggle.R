library(keras)

#-----------read csv file as compounds--------------
compounds <- read.csv('../input/traintakeda/train.csv')

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,-c(1)]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]

#-----------read csv file as compounds--------------
compounds.t <- read.csv('../input/testtakeda/test.csv')

#-----------remove some columns if needed--------------
trimed.compounds.t <- compounds.t[,-c(1)]

#-----------select rows without empty cells---------
is.completes.t <- complete.cases(trimed.compounds.t)
is.completes.t

complete.compounds.t <- trimed.compounds.t[is.completes.t,]

#--------------------divide into test and training data----------------------
train_size = 1

n = nrow(complete.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-----------select x from the dataset-----------------
x.0 <- complete.compounds[,-c(1)]
x.0.t <- complete.compounds.t[,]

x.sds <- apply(x.0, 2, sd)
x.sds.t <- apply(x.0.t, 2, sd)

sd.is.not.0 <- x.sds != 0 
sd.is.not.0
sd.is.not.0.t <- x.sds.t != 0
sd.is.not.0.t
x.0 <- x.0[, (sd.is.not.0) * (sd.is.not.0.t)]
x.0.t <- x.0.t[, (sd.is.not.0) * (sd.is.not.0.t)]

#pc = 100
#rpca = prcomp(x = x.0,scale=T)
#x.0 <- as.data.frame(rpca$x[,c(1:pc)])

#rpca.t = prcomp(x = x.0.t,scale=T)
#x.0.t <- as.data.frame(rpca.t$x[,c(1:pc)])

#sampledata <- as.data.frame(cbind(train_labels, train_data))
#lsample <- lm(train_labels~., data = sampledata)
#summary(lsample)
#p <- predict(lsample)
#plot(p, train_labels)
#train_labels
train_data <- cbind.data.frame(x.0[perm,])
test_data <- cbind.data.frame(x.0[-perm,])
#-----------select y from the dataset------------------
train_labels <- complete.compounds[perm,c(1)]
test_labels <- complete.compounds[-perm,c(1)]

# Test data is *not* used when calculating the mean and std.

# Normalize training data
train_data <- scale(train_data) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(train_data, "scaled:center") 
col_stddevs_train <- attr(train_data, "scaled:scale")
test_data <- scale(test_data, center = col_means_train, scale = col_stddevs_train)

train_data[1, ] # First training sample, normalized

#-----------select x from the dataset-----------------
test_data.t <- cbind.data.frame(x.0.t[,])
# Test data is *not* used when calculating the mean and std.
# Use means and standard deviations from training set to normalize test set
test_data.t <- scale(test_data.t, center = col_means_train, scale = col_stddevs_train)

#--------------modeling------------------------
build_model <- function() {
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 300, activation = "relu",
                input_shape = dim(train_data)[2]) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.25) %>%
#    layer_dense(units = 100, activation = "relu") %>%
#    layer_batch_normalization() %>%
#    layer_dropout(rate = 0.25) %>%
    layer_dense(units = 1)
  
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adadelta(lr = initial_lr),
    metrics = list("mean_absolute_error")
  )
  
  model
}

epochs <- 300
batch_size <- 128
initial_lr<-0.01
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
  callbacks = list(print_dot_callback, cb_lr)
)

plot(history)

train_predictions <- model %>% predict(train_data)
plot(train_predictions, train_labels)
r2train <- cor(train_predictions, train_labels)**2
r2train

#-----------read csv file as compounds--------------
compounds.t <- read.csv('../input/testtakeda/test.csv')

#-----------remove some columns if needed--------------
trimed.compounds.t <- compounds.t[,-c(1)]

#-----------select rows without empty cells---------
is.completes.t <- complete.cases(trimed.compounds.t)
is.completes.t

complete.compounds.t <- trimed.compounds.t[is.completes.t,]

#-----------select x from the dataset-----------------
x.0.t <- complete.compounds.t[,]
x.0.t <- x.0.t[, sd.is.not.0]

rpca = prcomp(x = x.0.t,scale=T)
x.0.t <- as.data.frame(rpca$x[,c(1:pc)])

test_data.t <- cbind.data.frame(x.0.t[,])
# Test data is *not* used when calculating the mean and std.
# Use means and standard deviations from training set to normalize test set
test_data.t <- scale(test_data.t, center = col_means_train, scale = col_stddevs_train)

test_predictions.t <- model %>% predict(test_data.t)
