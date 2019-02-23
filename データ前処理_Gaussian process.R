#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)
View(compounds)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]
View(complete.compounds)

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(11:14)]
View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)
x.sds[1]

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
View(x)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]
y

#-----------standarization of y------------------------
y <- (y - mean(y)) / sd(y)
preprocessed.y <- (y - min(y)) / (max(y) - min(y))
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x------------------------
x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
preprocessed.x <- apply(x, 2, function(x) {(x - min(x)) / (max(x) - min(x))})

View(preprocessed.x)

#-----------x converted into data frame type for machine learning-----------
#class(preprocessed.x)
preprocessed.x0 <- data.frame(preprocessed.x[,])
#class(preprocessed.x)

#apply(preprocessed.x, 2, mean)
#apply(preprocessed.x, 2, sd)

preprocessed.x1 <- preprocessed.x0[c(1:100),]
preprocessed.x2 <- preprocessed.x0[c(101:200),]
preprocessed.y1 <- preprocessed.y[1:100]
preprocessed.y2 <- preprocessed.y[101:200]
