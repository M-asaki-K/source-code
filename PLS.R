#-----------definition of scaled.compounds (used for PLS)-----------
scaled.compounds <- cbind(preprocessed.y, preprocessed.x)
View(scaled.compounds)
scaled.compounds <- as.data.frame(scaled.compounds)

#-------ここまでデータ前処理--------
#PLS
library(pls)
help(package="pls")

compounds.plsr <- plsr(preprocessed.y~., data=scaled.compounds, validation="CV")
summary(compounds.plsr)
plot(compounds.plsr, "validation")

compounds.plsr$coefficients[, , 5]
predict(compounds.plsr)
predict(compounds.plsr)[, , 5]
plsr.predicted.y <- predict(compounds.plsr)[, , 5]
plsr.r2 <- cor(preprocessed.y, plsr.predicted.y)**2
plsr.r2

plot(preprocessed.y, plsr.predicted.y,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLSR")
abline(a=0, b=1)
