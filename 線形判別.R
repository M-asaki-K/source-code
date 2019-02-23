nominal.y <- cut(compounds[, 2], breaks=c(-Inf, 1500, Inf), labels=c(1, 2))
nominal.y

library(MASS)
help(lda)

table(nominal.y)

category.compounds <- cbind(nominal.y, multi.regression.x)
category.compounds<-as.data.frame(category.compounds)
compounds.lda <- lda(nominal.y~., data=category.compounds)
compounds.lda

lda.predicted.y <- predict(compounds.lda)$class

plot(compounds.lda)

lda.table <- table(nominal.y, lda.predicted.y)
lda.table

sample.size <- sum(lda.table)
sample.size

lda.tp <- lda.table[2, 2]
lda.tp

lda.tn <- lda.table[1, 1]
lda.tn

lda.correct.ratio <- (lda.tp + lda.tn) / sample.size * 100
lda.correct.ratio

positive.sample.size <- sum(lda.table[2, ])
positive.sample.size

lda.sensitivity <- lda.tp / positive.sample.size * 100
lda.sensitivity

