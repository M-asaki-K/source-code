#ˆÈ‰ºstepwise‚É‚æ‚éAICÅ¬‰»
compounds.lm <- lm(preprocessed.y~., data=multi.regression.compounds)
compounds.lm
step(compounds.lm)
