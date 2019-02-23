library(gaselect)

ctrl <- genAlgControl(populationSize = 100, numGenerations = 30, minVariables = 5,
                      maxVariables = 18, verbosity = 1)

evaluatorSRCV <- evaluatorPLS(numReplications = 2, innerSegments = 7, testSetSize = 0.4,
                              numThreads = 1)

evaluatorRDCV <- evaluatorPLS(numReplications = 2, innerSegments = 5, outerSegments = 3,
                              numThreads = 1)

# Generate demo-data
set.seed(12345)
X <- as.matrix(preprocessed.x)
#View(X)
y <- (preprocessed.y)
y

resultSRCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorSRCV, seed = 123)
resultRDCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorRDCV, seed = 123)

subsets(resultSRCV, 1:5)
subsets(resultRDCV, 1:5)


