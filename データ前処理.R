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
x <- complete.compounds[,c(3:35)]
View(x)
x.2 <- complete.compounds[,c(6:11,13,17:25)]
View(x.2)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)
x.sds[1]

x.mean <- apply(x, 2, mean)
x.mean[2]
x.mean

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
View(x)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(2)]
y

#-----------standarization of y------------------------
preprocessed.y <- (y - mean(y)) / sd(y)
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x------------------------
apply(x, 2, mean)
apply(x, 2, sd)
preprocessed.x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
View(preprocessed.x)

#-----------x converted into data frame type for machine learning-----------
#class(preprocessed.x)
preprocessed.x <- data.frame(preprocessed.x)
#class(preprocessed.x)

#apply(preprocessed.x, 2, mean)
#apply(preprocessed.x, 2, sd)
