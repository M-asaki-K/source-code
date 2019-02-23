#-----------compare the number of columns and rows--------------
ncol(preprocessed.x)
nrow(preprocessed.x)

#-----------pick up columns if needed---------------------------
multi.regression.x <- preprocessed.x[ , ]
View(multi.regression.x)

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(preprocessed.y, multi.regression.x)
View(multi.regression.compounds)

#-----------verify each factor's relations--------------------------
pairs(multi.regression.compounds)
cor(multi.regression.compounds)
multi.regression.compounds <- as.data.frame(multi.regression.compounds)

#-----------MLR regression--------------------------------------------
compounds.lm <- lm(preprocessed.y~., data=multi.regression.compounds)
compounds.lm

summary(compounds.lm)

lm.predicted.y <- predict(compounds.lm)
lm.predicted.y

cor(preprocessed.y, lm.predicted.y)
lm.r2 <- cor(preprocessed.y, lm.predicted.y)**2
lm.r2

plot(preprocessed.y, lm.predicted.y)
abline(a=0, b=1)

plot(preprocessed.y, lm.predicted.y,
     xlab="Observed value",
     ylab="Predicted value")
abline(a=0, b=1)

