pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#if you want to change the number of threads for the calculation, please change the value "detectCores()"
registerDoParallel(makeCluster(detectCores()))
library(rcdk)
library(rcdklibs)

#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)
View(compounds)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,c(2:4)] #choose smiles, EA and IE

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]

#---------------Calculating molecules features------------------
params <- detectCores()
df2 <- complete.compounds
df2

### SPLIT DATA INTO K FOLDS ###
set.seed(2016)
df2$fold <- caret::createFolds(1:nrow(df2), k = params, list = FALSE)
max(df2$fold)

out <- foreach(j = 1:params, .combine = rbind, .inorder = FALSE, .packages = c("rcdk", "rcdklibs")) %dopar% {
  test <- df2[df2$fold == j, ]
  mols <- parse.smiles(as.character(test[, c(1)]))
  dnames <- get.desc.names('constitutional')
  dnames2 <- get.desc.names('electronic')
  dnames3 <- get.desc.names('hybrid')
  dnames4 <- get.desc.names('topological')
  dnames5 <- get.desc.names('geometrical')
  descs <- eval.desc(mols, cbind(dnames,dnames2,dnames3,dnames4,dnames5), verbose=FALSE)
  data.frame(descs)
}

tx<-t(out)
txomit<-na.omit(tx) #omitting columns of NAs
Alldescriptors<-t(txomit)

x <- as.data.frame(Alldescriptors)
View(x)
