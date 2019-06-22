devtools::install_github("rstudio/keras")
devtools::install_github("rstudio/tensorflow")

library(keras)
install_keras()
install_tensorflow
mnist = dataset_mnist()

install.packages('keras', dependencies=T)
install_keras()

#学習用説明変数の作成
x_train = mnist$train$x
#学習用目的変数の作成
y_train = mnist$train$y
#テスト用説明変数の作成
x_test = mnist$test$x
#テスト用目的変数の作成
y_test = mnist$test$y

mnist$train$x[1,,]

#学習用目的変数の頭出し
head(mnist$train$y)
#学習用目的変数の構成要素
unique(mnist$train$y)

#元の次元数
dim(x_train)
dim(x_test)
#各データ28行28列の行列のところを
#各データ1行784列に変換します
dim(x_train) = c(nrow(x_train), 784)
dim(x_test) = c(nrow(x_test), 784)
#変換後の次元を確認します。
dim(x_train)
dim(x_test)

#x_trainの最大値
max(x_train)
#x_trainの最小値
min(x_train)
#x_testの最大値
max(x_test)
#x_testの最小値
min(x_test)

#説明変数の基準化
x_train = x_train / 255
x_test = x_test / 255

#基準化の結果
#x_trainの最大値
max(x_train)
#x_trainの最小値
min(x_train)
#x_testの最大値
max(x_test)
#x_testの最小値
min(x_test)

#目的変数のone-hotベクトル化
y_train = to_categorical(y_train, 10)
y_test = to_categorical(y_test, 10)
#目的変数変換の結果(それぞれ一つ目の成分)
y_train[1,]
y_test[1,]

#Sequential Modelの呼び出し
model <- keras_model_sequential()
#784→32→relu→10→softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dense(units = 10) %>%
  layer_activation('softmax')

summary(model)

#モデルのコンパイル
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#学習と評価
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#テストデータを使った最終的な評価の結果
model %>% evaluate(x_test, y_test)

#テストデータを使った予測値の算出
prediction=model %>% predict_classes(x_test)
head(prediction)

#Sequential Modelの呼び出し
model = keras_model_sequential()
#784→32→relu→Dropout(0.3)→10→Dropout(0.3)→softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10) %>%
  layer_dropout(rate = 0.3) %>%
  layer_activation('softmax')

#summary関数でNeuralNetwork構造確認
summary(model)

#モデルのコンパイル
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#学習と評価
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#lossとaccのプロット
plot(history)

#テストデータを使った最終的な評価の結果
model %>% evaluate(x_test, y_test)

#Sequential Modelの呼び出し
model = keras_model_sequential()
#784→32→relu→Dropout(0.3)→10→Dropout(0.3)→softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10) %>%
  layer_dropout(rate = 0.1) %>%
  layer_activation('softmax')

#summary関数でNeuralNetwork構造確認
summary(model)

#モデルのコンパイル
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#学習と評価
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#lossとaccのプロット
plot(history)

#テストデータを使った最終的な評価の結果
model %>% evaluate(x_test, y_test)

#Sequential Modelの呼び出し
model = keras_model_sequential()

#784→256→relu→Dropout(0.4)→128→relu→Dropout(0.3)→10→softmax
model %>%
  layer_dense(units = 256, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 128) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10) %>%
  layer_activation('softmax')

#summary関数でNeuralNetwork構造確認
summary(model)

#モデルのコンパイル
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#学習と評価
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#lossとaccのプロット
plot(history)

#テストデータを使った最終的な評価の結果
model %>% evaluate(x_test, y_test)

batch_size = 128
num_classes = 10
epochs = 5
#epochsは5くらいの方が早く終わってよいと思います。
#epochs12だと少なくとも60分かかります。

#入力次元
img_rows = 28
img_cols = 28

#データのロードと学習用データとテスト用データの作成
mnist = dataset_mnist()
x_train = mnist$train$x
y_train = mnist$train$y
x_test = mnist$test$x
y_test = mnist$test$y

dim(x_train) = c(nrow(x_train), img_rows, img_cols, 1)
dim(x_test) = c(nrow(x_test), img_rows, img_cols, 1)
input_shape = c(img_rows, img_cols, 1)

x_train = x_train / 255
x_test = x_test / 255

#目的変数をone-hotベクトルへ変換
y_train = to_categorical(y_train, num_classes)
y_test = to_categorical(y_test, num_classes)

#モデルの定義
model = keras_model_sequential()
model %>%
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = input_shape) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_dropout(rate = 0.25) %>%
  layer_flatten() %>%
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = num_classes, activation = 'softmax')

#モデルのコンパイル
model %>% compile(
  loss = loss_categorical_crossentropy,
  optimizer = optimizer_adadelta(),
  metrics = c('accuracy')
)

#学習と評価
model %>% fit(
  x_train, y_train,
  batch_size = batch_size,
  epochs = epochs,
  verbose = 1,
  validation_data = list(x_test, y_test)
)
#テストデータでの評価
model %>% evaluate(x_test, y_test)
