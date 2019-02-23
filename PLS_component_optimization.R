#PLS
library(pls)
help(package="pls")

compounds.plsr <- plsr(preprocessed.y~., data=scaled.compounds, validation="CV")
ncomp.onesigma <- selectNcomp(compounds.plsr, method = "onesigma", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma

predict(compounds.plsr)[, , ncomp.onesigma]
plsr.predicted.y <- predict(compounds.plsr)[, , ncomp.onesigma]
plsr.r2 <- cor(preprocessed.y, plsr.predicted.y)**2
plsr.r2

plot(preprocessed.y, plsr.predicted.y,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLSR")
abline(a=0, b=1)

compounds.plsr$coefficients[, , ncomp.onesigma]
