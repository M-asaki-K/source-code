#---------------------get SMILES and evaluate descriptors-------------------
smiles <- c('CCC', 'c1ccccc1', 'CC(=O)C','CC[Al](CC)CC')
mols <- sapply(smiles, parse.smiles)

dnames <- get.desc.names('constitutional')
dnames2 <- get.desc.names('electronic')
dnames3 <- get.desc.names('hybrid')
dnames4 <- get.desc.names('topological')
dnames5 <- get.desc.names('geometrical')
descs <- eval.desc(mols, cbind(dnames,dnames2,dnames3,dnames4,dnames5), verbose=TRUE)

View(descs)

#------------------------bp regression--------------------------
data(bpdata)
mols <- parse.smiles(bpdata[,1])
#------------------descriptors manually selected----------------------
descNames <- c(
  'org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor',
  'org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor',
  'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor')
descs <- eval.desc(mols, descNames)

#-----------------manual selection of x------------------------
model <- lm(BP ~ khs.sCH3 + khs.sF + apol + nHBDon, data.frame(bpdata, descs))
summary(model)

pr <- predict(model,descs)
plot(bpdata[,c(2)],pr)
