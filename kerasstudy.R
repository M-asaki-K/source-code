devtools::install_github("rstudio/keras")
devtools::install_github("rstudio/tensorflow")

library(keras)
install_keras()
install_tensorflow
mnist = dataset_mnist()

install.packages('keras', dependencies=T)
install_keras()

#�w�K�p�����ϐ��̍쐬
x_train = mnist$train$x
#�w�K�p�ړI�ϐ��̍쐬
y_train = mnist$train$y
#�e�X�g�p�����ϐ��̍쐬
x_test = mnist$test$x
#�e�X�g�p�ړI�ϐ��̍쐬
y_test = mnist$test$y

mnist$train$x[1,,]

#�w�K�p�ړI�ϐ��̓��o��
head(mnist$train$y)
#�w�K�p�ړI�ϐ��̍\���v�f
unique(mnist$train$y)

#���̎�����
dim(x_train)
dim(x_test)
#�e�f�[�^28�s28��̍s��̂Ƃ����
#�e�f�[�^1�s784��ɕϊ����܂�
dim(x_train) = c(nrow(x_train), 784)
dim(x_test) = c(nrow(x_test), 784)
#�ϊ���̎������m�F���܂��B
dim(x_train)
dim(x_test)

#x_train�̍ő�l
max(x_train)
#x_train�̍ŏ��l
min(x_train)
#x_test�̍ő�l
max(x_test)
#x_test�̍ŏ��l
min(x_test)

#�����ϐ��̊��
x_train = x_train / 255
x_test = x_test / 255

#����̌���
#x_train�̍ő�l
max(x_train)
#x_train�̍ŏ��l
min(x_train)
#x_test�̍ő�l
max(x_test)
#x_test�̍ŏ��l
min(x_test)

#�ړI�ϐ���one-hot�x�N�g����
y_train = to_categorical(y_train, 10)
y_test = to_categorical(y_test, 10)
#�ړI�ϐ��ϊ��̌���(���ꂼ���ڂ̐���)
y_train[1,]
y_test[1,]

#Sequential Model�̌Ăяo��
model <- keras_model_sequential()
#784��32��relu��10��softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dense(units = 10) %>%
  layer_activation('softmax')

summary(model)

#���f���̃R���p�C��
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#�w�K�ƕ]��
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#�e�X�g�f�[�^���g�����ŏI�I�ȕ]���̌���
model %>% evaluate(x_test, y_test)

#�e�X�g�f�[�^���g�����\���l�̎Z�o
prediction=model %>% predict_classes(x_test)
head(prediction)

#Sequential Model�̌Ăяo��
model = keras_model_sequential()
#784��32��relu��Dropout(0.3)��10��Dropout(0.3)��softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10) %>%
  layer_dropout(rate = 0.3) %>%
  layer_activation('softmax')

#summary�֐���NeuralNetwork�\���m�F
summary(model)

#���f���̃R���p�C��
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#�w�K�ƕ]��
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#loss��acc�̃v���b�g
plot(history)

#�e�X�g�f�[�^���g�����ŏI�I�ȕ]���̌���
model %>% evaluate(x_test, y_test)

#Sequential Model�̌Ăяo��
model = keras_model_sequential()
#784��32��relu��Dropout(0.3)��10��Dropout(0.3)��softmax
model %>%
  layer_dense(units = 32, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 10) %>%
  layer_dropout(rate = 0.1) %>%
  layer_activation('softmax')

#summary�֐���NeuralNetwork�\���m�F
summary(model)

#���f���̃R���p�C��
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#�w�K�ƕ]��
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#loss��acc�̃v���b�g
plot(history)

#�e�X�g�f�[�^���g�����ŏI�I�ȕ]���̌���
model %>% evaluate(x_test, y_test)

#Sequential Model�̌Ăяo��
model = keras_model_sequential()

#784��256��relu��Dropout(0.4)��128��relu��Dropout(0.3)��10��softmax
model %>%
  layer_dense(units = 256, input_shape = c(784)) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 128) %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10) %>%
  layer_activation('softmax')

#summary�֐���NeuralNetwork�\���m�F
summary(model)

#���f���̃R���p�C��
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

#�w�K�ƕ]��
history = model %>% fit(
  x_train, y_train,
  epochs = 30, batch_size = 128,
  validation_split = 0.2
)

#loss��acc�̃v���b�g
plot(history)

#�e�X�g�f�[�^���g�����ŏI�I�ȕ]���̌���
model %>% evaluate(x_test, y_test)

batch_size = 128
num_classes = 10
epochs = 5
#epochs��5���炢�̕��������I����Ă悢�Ǝv���܂��B
#epochs12���Ə��Ȃ��Ƃ�60��������܂��B

#���͎���
img_rows = 28
img_cols = 28

#�f�[�^�̃��[�h�Ɗw�K�p�f�[�^�ƃe�X�g�p�f�[�^�̍쐬
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

#�ړI�ϐ���one-hot�x�N�g���֕ϊ�
y_train = to_categorical(y_train, num_classes)
y_test = to_categorical(y_test, num_classes)

#���f���̒�`
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

#���f���̃R���p�C��
model %>% compile(
  loss = loss_categorical_crossentropy,
  optimizer = optimizer_adadelta(),
  metrics = c('accuracy')
)

#�w�K�ƕ]��
model %>% fit(
  x_train, y_train,
  batch_size = batch_size,
  epochs = epochs,
  verbose = 1,
  validation_data = list(x_test, y_test)
)
#�e�X�g�f�[�^�ł̕]��
model %>% evaluate(x_test, y_test)